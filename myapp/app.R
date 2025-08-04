# Shinylive fixes -------------------------------------------------------------------------------------------------

# 1. Build-time safety switches
options(
  bslib.download_fonts = FALSE, # never bundle Google-font files
  bslib.prefer_remote_fonts = TRUE # write plain @import rules if needed
)

# 2. Helper for wasm packages
if (!nzchar(Sys.getenv("SHINYLIVE"))) {
  install_missing_pkgs <- function(pkgs) {}
} else {
  install_missing_pkgs <- function(pkgs) {
    if (!requireNamespace("webr", quietly = TRUE)) {
      return(invisible())
    }
    need <- setdiff(pkgs, rownames(installed.packages()))
    if (length(need)) webr::install(need)
  }
}
install_missing_pkgs(c("rsimsum", "ggplot2"))


# Packages --------------------------------------------------------------------------------------------------------

library(shiny)
library(bslib)
library(DT)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(readr)
library(rsimsum)
library(munsell)

# Import data -----------------------------------------------------------------------------------------------------

ds_cors <- read_csv("data/cors_small.csv.gz", show_col_types = FALSE)
ds_energy <- read_csv("data/energy_small.csv.gz", show_col_types = FALSE)
ds_models <- read_csv("data/models_small.csv.gz", show_col_types = FALSE)

collect_df <- identity

# Helpers & plotting ----------------------------------------------------------------------------------------------

app_theme <- bs_theme(version = 5, preset = "bootstrap")

pretty_dgm <- function(x) recode(x, nb_bin = "NegBin", pois_bin = "Poisson")

sim_params <- ds_models %>%
  distinct(ate_exp, mech, phi, dgm, method, weight) %>%
  arrange(ate_exp, mech, phi)

format_rsimsum_summary <- function(simsum_obj) {
  if (is.null(simsum_obj)) {
    return(NULL)
  }
  target_stats <- c("bias", "rbias", "empse", "modelse", "cover")
  grouping_vars <- if (is.null(simsum_obj$by)) character(0) else simsum_obj$by
  
  wide <- simsum_obj$summ %>%
    filter(stat %in% target_stats) %>%
    group_by(across(all_of(c("method", grouping_vars, "stat")))) %>%
    summarise(est = mean(est), mcse = mean(mcse), .groups = "drop") %>%
    pivot_wider(
      id_cols = c(all_of(grouping_vars), method),
      names_from = stat, values_from = c(est, mcse),
      names_glue = "{stat}_{.value}"
    )
  
  exp_cols <- as.vector(outer(target_stats, c("est", "mcse"), paste, sep = "_"))
  wide[setdiff(exp_cols, names(wide))] <- NA_real_
  
  if ("dgm" %in% grouping_vars) wide$dgm <- pretty_dgm(wide$dgm)
  
  wide %>%
    mutate(
      `Bias (MCSE)` = sprintf("%.3f (%.3f)", bias_est, bias_mcse),
      `Rel. Bias (MCSE)` = sprintf("%.1f%% (%.1f%%)", rbias_est * 100, rbias_mcse * 100),
      `Emp. SE` = sprintf("%.3f", empse_est),
      `Model SE` = sprintf("%.3f", modelse_est),
      `Coverage (MCSE)` = sprintf("%.3f (%.3f)", cover_est, cover_mcse)
    ) %>%
    select(all_of(grouping_vars),
           Method = method,
           `Bias (MCSE)`, `Rel. Bias (MCSE)`, `Emp. SE`, `Model SE`, `Coverage (MCSE)`
    ) |>
    rename('Exposure' = dgm)
}

plot_zip_with_estimates <- function(data, true_value) {
  data$model <- data$dgm
  data$dgm <- pretty_dgm(data$dgm)
  
  data <- data %>%
    filter(!is.na(std.error) & std.error > 0) %>%
    mutate(
      z_score = (estimate - true_value) / std.error,
      abs_z_score = abs(z_score),
      includes_true_value =
        ifelse(conf.low <= true_value & conf.high >= true_value, "Yes", "No")
    ) %>%
    group_by(method) %>%
    mutate(centile = percent_rank(abs_z_score)) %>%
    ungroup()
  
  if (nrow(data) == 0) {
    return(ggplot() +
             labs(title = "No data.") +
             theme_void())
  }
  
  perc95 <- data %>%
    group_by(method) %>%
    summarise(p95 = quantile(centile, .95), .groups = "drop")
  
  ggplot(data, aes(x = estimate, y = centile)) +
    geom_errorbarh(aes(
      xmin = conf.low, xmax = conf.high,
      colour = includes_true_value
    ), height = 0.002) +
    geom_point(colour = "white", size = 0.7) +
    geom_vline(xintercept = true_value, colour = "red", linetype = "dashed") +
    scale_color_manual(values = c(Yes = "black", No = "orange")) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    labs(
      x = "Estimate", y = "Fractional centile of |Z-score|",
      colour = "CI includes\ntrue value"
    ) +
    theme_minimal() +
    facet_grid(dgm ~ method) +
    geom_hline(
      data = perc95, aes(yintercept = p95),
      colour = "blue", linetype = "dotted"
    ) +
    ylim(0.5, 1)
}

balance_summary_tbl <- function(df_energy, df_cors) {
  ess_tbl <- df_energy %>%
    rename(eps_A = energy_A, eps_C = energy_X) %>%
    group_by(dgm, method) %>%
    summarise(
      Mean_ESS = mean(ess), SD_ESS = sd(ess),
      P5_ESS = quantile(ess, .05), P95_ESS = quantile(ess, .95),
      Dw_mean = mean(D_w), Dw_sd = sd(D_w),
      eps_A_mean = mean(eps_A), eps_C_mean = mean(eps_C),
      .groups = "drop"
    )
  
  cor_tbl <- df_cors %>%
    mutate(abs_cor = abs(cor)) %>%
    group_by(dgm, method, rep) %>%
    summarise(mean_abs_cor = mean(abs_cor), .groups = "drop") %>%
    group_by(dgm, method) %>%
    summarise(
      Mean_rho = mean(mean_abs_cor), SD_rho = sd(mean_abs_cor),
      P95_rho = quantile(mean_abs_cor, .95), Max_rho = max(mean_abs_cor),
      .groups = "drop"
    )
  
  ess_tbl %>%
    inner_join(cor_tbl, by = c("dgm", "method")) %>%
    select(dgm,
           Method = method,
           Mean_ESS, SD_ESS, P5_ESS, P95_ESS,
           Dw_mean, Dw_sd, eps_A_mean, eps_C_mean,
           Mean_rho, SD_rho, P95_rho, Max_rho
    ) %>%
    arrange(dgm, Method) |>
    rename(Exposure = dgm)
}

round_numeric_cols <- function(df) {
  num_cols <- names(df)[vapply(df, is.numeric, logical(1))]
  for (cl in num_cols) {
    x <- df[[cl]]
    ch <- as.character(x)
    tiny <- abs(x) < 1e-4 & x != 0
    small <- abs(x) < 1e-3 & abs(x) >= 1e-4
    ch[tiny] <- "<0.0001"
    ch[small] <- "<0.001"
    keep <- !(tiny | small)
    if (any(keep)) {
      if (grepl("ESS", cl)) {
        ch[keep] <- format(round(x[keep]), scientific = FALSE, trim = TRUE)
      } else {
        rng <- range(x[keep])
        ch[keep] <- if (max(abs(rng)) < 1e3 && max(abs(rng)) > 1e-3) {
          sprintf("%.3f", round(x[keep], 3))
        } else {
          format(signif(x[keep], 3), scientific = TRUE)
        }
      }
    }
    df[[cl]] <- ch
  }
  df
}

prettify_headers <- function(df) {
  nm <- names(df)
  map <- c(
    Mean_ESS = "Mean ESS",
    SD_ESS = "SD ESS",
    P5_ESS = "5<sup>th</sup> perc ESS",
    P95_ESS = "95<sup>th</sup> perc ESS",
    Dw_mean = "Dw",
    Dw_sd = "SD Dw",
    eps_A_mean = "&epsilon;<sub>A</sub>",
    eps_C_mean = "&epsilon;<sub>C</sub>",
    Mean_rho = paste0(
      "<span style='display:inline-block;",
      " padding:0 0.4em;",               # space for the caret
      " border-left:1px solid currentColor;",
      " border-right:1px solid currentColor;'>",
      "&rho;<sub>w</sub></span>"
    ),
    SD_rho = paste0(
      "<span style='display:inline-block;",
      " padding:0 0.4em;",
      " border-left:1px solid currentColor;",
      " border-right:1px solid currentColor;'>",
      "&rho;<sub>w</sub></span>"
    ),
    P95_rho = paste0(
      "95<sup>th</sup> ",
      "<span style='display:inline-block;",
      " padding:0 0.4em;",
      " border-left:1px solid currentColor;",
      " border-right:1px solid currentColor;'>",
      "&rho;<sub>w</sub></span>"
    ),
    Max_rho = paste0(
      "Max ",
      "<span style='display:inline-block;",
      " padding:0 0.4em;",
      " border-left:1px solid currentColor;",
      " border-right:1px solid currentColor;'>",
      "&rho;<sub>w</sub></span>"
    )
  )
  names(df) <- ifelse(nm %in% names(map), map[nm], nm)
  df
}

source("plot_R/coverage_plot.R")  ### NEW: source external plotting helper

ui <- fluidPage(
  ## Roboto at run-time (no build-time effect) -------------------------------
  tags$head(
    tags$link(
      rel  = "stylesheet",
      href = "https://fonts.googleapis.com/css2?family=Roboto:wght@400;700&display=swap"
    ),
    tags$style(HTML(":root{--bs-body-font-family:'Roboto',system-ui,sans-serif;}")),
    
    ## jsPDF from CDN for PDF creation
    tags$script(src = "https://cdnjs.cloudflare.com/ajax/libs/jspdf/2.5.1/jspdf.umd.min.js"),
    
    ## Custom message handlers for CSV & PDF downloads
    tags$script(HTML("
      Shiny.addCustomMessageHandler('download-csv', function(msg) {
        var blob = new Blob([msg.data], { type: 'text/csv;charset=utf-8;' });
        var a = document.createElement('a');
        var url = URL.createObjectURL(blob);
        a.href = url; a.download = msg.filename;
        document.body.appendChild(a); a.click();
        document.body.removeChild(a); URL.revokeObjectURL(url);
      });
      Shiny.addCustomMessageHandler('download-pdf', function(msg) {
        const { jsPDF } = window.jspdf;
        // Grab the <img> for the currently visible plot (Zip or Coverage):
        function findVisiblePlotImage(ids) {
          for (const id of ids) {
            var container = document.getElementById(id);
            if (!container) continue;
            var img = container.getElementsByTagName('img')[0];
            if (!img) continue;
            if (container.offsetParent !== null &&
                window.getComputedStyle(container).display !== 'none' &&
                window.getComputedStyle(container).visibility !== 'hidden') {
              return img;
            }
          }
          return null;
        }
        var img = findVisiblePlotImage(['zip_plot','coverage_plot']) ||
                  (document.getElementById('zip_plot')||{}).getElementsByTagName?.('img')[0] ||
                  (document.getElementById('coverage_plot')||{}).getElementsByTagName?.('img')[0];

        if (!img) { alert('Plot image not found'); return; }
        var doc = new jsPDF();
        // full-width image, preserve aspect ratio
        doc.addImage(img.src, 'PNG', 10, 10, 190, 0);
        doc.save(msg.filename);
      });
    "))
  ),
  
  theme = app_theme,
  
  tags$style(HTML("
  div.dataTables_filter{ float:left!important; text-align:left!important; }
  table.dataTable tbody td{ padding:6px 10px; }

  /* NEW ─ header cells: centred text & extra space */
  table.dataTable thead th{
    text-align:center;         /* centre content */
    vertical-align:bottom;     /* keep maths symbols on the baseline */
    padding:6px 16px 4px 16px; /* ↑ top | → right | ↓ bottom | ← left */
    white-space:nowrap;        /* stop long labels wrapping */
  }
")),
  
  titlePanel("Interactive Simulation Summaries"),
  
  sidebarLayout(
    sidebarPanel(
      width = 3,
      
      h4("Simulation Parameters"),
      selectInput("ate_exp", "Effect Size (ATE):", choices = NULL),
      selectInput("mech", "Missingness Mechanism:", choices = NULL),
      uiOutput("phi_ui"),
      selectInput("weight_type", "Winsorise at 99th perc:", 
                  choices = c("No" = "raw", "Yes" = "win"), selected = "raw"),
      checkboxGroupInput("dgm",  "Exposure Distribution:",       choices = NULL),
      
      hr(),
      h4("View Selection"),
      selectInput("view_type", "Result Type to View:",
                  choices = c("Model Performance" = "models",
                              "Covariate Balance"   = "balance")),
      
      hr(),
      uiOutput("method_ui"),
      
      hr(),
      h4("Downloads"),
      actionButton("btn_dl_csv", "Download table (.csv)",
                   style = "width:100%; margin-bottom:5px;"),
      actionButton("btn_dl_pdf", "Download plot (.pdf)",
                   style = "width:100%;")
    ),
    
    mainPanel(
      width = 9,
      tabsetPanel(
        id = "main_tabs",
        tabPanel("Summary Table", DTOutput("summary_tbl")),
        tabPanel("Zip Plot",      plotOutput("zip_plot", height = "800px")),
        tabPanel("Coverage Plot", plotOutput("coverage_plot", height = "800px"))  ### NEW
      )
    )
  )
)


server <- function(input, output, session) {
  ## Populate selectors -------------------------------------------------------
  updateSelectInput(session, "ate_exp",
                    choices  = setNames(sort(unique(sim_params$ate_exp)),
                                        paste0("ln ", sort(unique(sim_params$ate_exp)))),
                    selected = 1.1
  )
  updateSelectInput(session, "mech", choices = unique(sim_params$mech))
  updateCheckboxGroupInput(
    session, "dgm",
    choices  = setNames(unique(sim_params$dgm), pretty_dgm(unique(sim_params$dgm))),
    selected = "nb_bin"
  )
  
  output$phi_ui <- renderUI({
    req(input$mech)
    if (input$mech == "complete") {
      selectInput("phi", "Proportion of Missingness:", choices = 0, selected = 0)
    } else {
      selectInput("phi", "Proportion of Missingness:",
                  choices = sort(unique(sim_params$phi[sim_params$mech == input$mech])))
    }
  })
  
  output$method_ui <- renderUI({
    req(input$view_type)
    available <- if (input$view_type == "models") {
      unique(sim_params$method)
    } else {
      setdiff(unique(sim_params$method), c("unadjusted", "adjusted"))
    }
    checkboxGroupInput("method", "Method:", choices = available, selected = available)
  })
  
  ## Reactive filters ---------------------------------------------------------
  base_filter <- function(ds) {
    ds %>% filter(
      ate_exp == input$ate_exp,
      mech == input$mech,
      phi == req(input$phi),
      weight  == input$weight_type,
      dgm %in% input$dgm,
      method %in% input$method
    )
  }
  filt_energy <- reactive(collect_df(base_filter(ds_energy)) %>% mutate(dgm = pretty_dgm(dgm)))
  filt_cors   <- reactive(collect_df(base_filter(ds_cors))   %>% mutate(dgm = pretty_dgm(dgm)))
  filt_models <- reactive(collect_df(base_filter(ds_models)) %>% mutate(dgm = pretty_dgm(dgm)))
  
  ## Summary table ------------------------------------------------------------
  current_table_raw <- reactive({
    if (input$view_type == "models") {
      df <- filt_models()
      if (!nrow(df)) return(NULL)
      sims <- rsimsum::simsum(
        df, "estimate", "std.error",
        true = log(as.numeric(input$ate_exp)),
        methodvar = "method", by = "dgm"
      )
      format_rsimsum_summary(sims)
    } else {
      de <- filt_energy(); dc <- filt_cors()
      if (!nrow(de) || !nrow(dc)) return(NULL)
      balance_summary_tbl(de, dc)
    }
  })
  
  current_table_disp <- reactive({
    tbl <- current_table_raw()
    if (is.null(tbl)) return(NULL)
    tbl <- round_numeric_cols(tbl)
    if (input$view_type == "balance") tbl <- prettify_headers(tbl)
    tbl
  })
  
  output$summary_tbl <- renderDT({
    tbl <- current_table_disp()
    if (is.null(tbl)) {
      datatable(
        data.frame(Message = "No data for this combination."),
        class    = "compact",
        rownames = FALSE
      )
    } else {
      
      ## ---------- FOOTNOTE TEXT & CAPTION (new) --------------------------- ##
      footnote <- if (input$view_type == "models") {
        paste(
          "Point estimates and confidence intervals were derived from weighted",
          "Poisson regression models with a sandwich estimator for variance.",
          "Unweighted outcome regressions (adjusted and unadjusted) are shown",
          "for comparison. All IPTW methods used stabilised weights.<br />(np)CBPS –",
          "(Non-Parametric) Covariate Balancing Propensity Scores; DGM –",
          "Data-Generating Mechanism; Emp. SE – Empirical Standard Error; GBM –",
          "Generalised Boosted Models; MCSE – Monte Carlo Standard Error;",
          "Model SE – Model-Based Standard Error; Rel. Bias – Relative Bias.",
          collapse = " "
        )
      } else {
        paste(
          "The ε<sub>C</sub> metric and the absolute treatment–covariate",
          "correlations were averaged across the three confounders.<br />(np)CBPS –",
          "(Non-Parametric) Covariate Balancing Propensity Scores; DGM –",
          "Data-Generating Mechanism; Dw – distance metric (optimised by the",
          "energy-balancing approach); ε<sub>A</sub> – energy distance between the",
          "weighted and unweighted marginal distributions of the exposure;",
          "ε<sub>C</sub> – energy distance between the weighted and unweighted",
          "marginal distributions of the covariates; ESS – Effective Sample Size;",
          "GBM – Generalised Boosted Models; MCSE – Monte Carlo Standard Error;",
          "NegBin – Negative Binomial; ρ<sub>w</sub> – average weighted",
          "treatment-covariate correlation; perc – percentile; SD – standard",
          "deviation.",
          collapse = " "
        )
      }                                             ### ADDED ###
      ## -------------------------------------------------------------------- ##
      
      datatable(
        tbl,
        class     = "stripe hover compact order-column row-border",
        rownames  = FALSE,
        filter    = "top",
        extensions = c("FixedColumns"),
        options    = list(
          scrollX       = TRUE,                  # turn on horizontal scrolling
          pageLength    = 20,
          autoWidth     = TRUE,
          fixedColumns  = list(leftColumns = 1)  # ← pin the first column
        ),
        selection = "none",
        escape    = FALSE,
        caption   = tags$caption(                   ### ADDED ###
          style = "caption-side: bottom; text-align: left; font-size: 0.85em;",
          HTML(footnote)
        )
      )
    }
  })
  
  ## Zip plot -----------------------------------------------------------------
  current_plot <- reactive({
    req(input$view_type == "models")
    df <- filt_models()
    if (!nrow(df)) return(NULL)
    plot_zip_with_estimates(df, true_value = log(as.numeric(input$ate_exp)))
  })
  
  output$zip_plot <- renderPlot({
    g <- current_plot(); req(!is.null(g))
    print(g)
  })
  
  ## Coverage plot ------------------------------------------------------------
  current_coverage_plot <- reactive({
    req(input$view_type == "models")
    df <- collect_df(base_filter(ds_models))
    if (!nrow(df)) return(NULL)
    plot_coverage_by_method(df, true_value = log(as.numeric(input$ate_exp)))
  })
  
  output$coverage_plot <- renderPlot({
    g <- current_coverage_plot(); req(!is.null(g))
    print(g)
  })
  
  observe({
    if (input$view_type != "models") {
      hideTab("main_tabs", "Zip Plot")
      hideTab("main_tabs", "Coverage Plot")    ### NEW
    } else {
      showTab("main_tabs", "Zip Plot")
      showTab("main_tabs", "Coverage Plot")    ### NEW
    }
  })
  
  ## CLIENT-SIDE DOWNLOADS ----------------------------------------------------
  # 1) Table as CSV
  observeEvent(input$btn_dl_csv, {
    tbl <- current_table_disp()
    validate(need(!is.null(tbl), "No data to download"))
    csv_lines <- c(
      paste(colnames(tbl), collapse = ","),
      apply(tbl, 1, function(r) paste(r, collapse = ","))
    )
    csv_text <- paste(csv_lines, collapse = "\r\n")
    session$sendCustomMessage(
      "download-csv",
      list(
        data     = csv_text,
        filename = paste0("summary-", Sys.Date(), ".csv")
      )
    )
  })
  
  # 2) Plot as PDF
  observeEvent(input$btn_dl_pdf, {
    # JS side will now pick the visible plot (Zip or Coverage)
    session$sendCustomMessage(
      "download-pdf",
      list(
        filename = paste0("plot-", Sys.Date(), ".pdf")
      )
    )
  })
}

# Run app ---------------------------------------------------------------------------------------------------------

shinyApp(ui, server)
