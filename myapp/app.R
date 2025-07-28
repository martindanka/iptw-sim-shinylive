
# ─────────────────────────────────────────────────────────────────────────────
# 0) Hint to shinylive about extra wasm packages (scanned, never executed)
# ─────────────────────────────────────────────────────────────────────────────
if (FALSE) {
  library(rsimsum)  # needed in the browser
  # Note: We deliberately do NOT include arrow/officer/flextable here,
  #       since those are for local (server) use only.
}

# ─────────────────────────────────────────────────────────────────────────────
# 1) Tiny helper: install wasm builds in the browser just in case
# ─────────────────────────────────────────────────────────────────────────────
install_missing_pkgs <- function(pkgs) {
  if (!nzchar(Sys.getenv("SHINYLIVE"))) return(invisible())  # not in browser
  if (!requireNamespace("webr", quietly = TRUE)) return(invisible())
  need <- setdiff(pkgs, rownames(installed.packages()))
  if (length(need)) webr::install(need)
}
install_missing_pkgs(c("rsimsum"))

# ─────────────────────────────────────────────────────────────────────────────
# 2) Libraries
# ─────────────────────────────────────────────────────────────────────────────
library(shiny)
library(bslib)
library(DT)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(rsimsum)   # available in wasm (bundled) and locally

is_browser <- nzchar(Sys.getenv("SHINYLIVE"))
if (!is_browser) {
  # Local-only, not available in the browser
  suppressPackageStartupMessages({
    library(arrow)
    # officer/flextable are only used for Word export locally; keep optional
    if (requireNamespace("officer", quietly = TRUE) &&
        requireNamespace("flextable", quietly = TRUE)) {
      # loaded on demand inside download handler
    }
  })
}

# Helper to safely materialise a data source whether it's Arrow or data.frame
collect_df <- function(x) {
  tryCatch(arrow::collect(x), error = function(e) x)
}

# ─────────────────────────────────────────────────────────────────────────────
# 3) Data: Arrow when local; .rds when in the browser
# ─────────────────────────────────────────────────────────────────────────────

ds_cors   <- readRDS("myapp/data/cors_small.rds")
ds_energy <- readRDS("myapp/data/energy_small.rds")
ds_models <- readRDS("myapp/data/models_small.rds")

# ─────────────────────────────────────────────────────────────────────────────
# 4) Helpers & plotting (adapted from your original code)
# ─────────────────────────────────────────────────────────────────────────────
pretty_dgm <- function(x) dplyr::recode(x, nb_bin = "NegBin", pois_bin = "Poisson")

sim_params <- collect_df(
  ds_models %>%
    distinct(ate_exp, mech, phi, dgm, method, weight) %>%
    arrange(ate_exp, mech, phi)
)

format_rsimsum_summary <- function(simsum_obj) {
  if (is.null(simsum_obj)) return(NULL)
  target_stats  <- c("bias", "rbias", "empse", "modelse", "cover")
  grouping_vars <- if (is.null(simsum_obj$by)) character(0) else simsum_obj$by
  
  wide <- simsum_obj$summ %>%
    filter(stat %in% target_stats) %>%
    group_by(across(all_of(c("method", grouping_vars, "stat")))) %>%
    summarise(est = mean(est), mcse = mean(mcse), .groups = "drop") %>%
    pivot_wider(id_cols = c(all_of(grouping_vars), method),
                names_from = stat, values_from = c(est, mcse),
                names_glue = "{stat}_{.value}")
  
  exp_cols <- as.vector(outer(target_stats, c("est", "mcse"), paste, sep = "_"))
  wide[setdiff(exp_cols, names(wide))] <- NA_real_
  
  if ("dgm" %in% grouping_vars) wide$dgm <- pretty_dgm(wide$dgm)
  
  wide %>%
    mutate(
      Bias        = sprintf("%.3f (%.3f)", bias_est,   bias_mcse),
      `Rel. Bias` = sprintf("%.1f%% (%.1f%%)", rbias_est * 100, rbias_mcse * 100),
      `Emp. SE`   = sprintf("%.3f", empse_est),
      `Model SE`  = sprintf("%.3f", modelse_est),
      Coverage    = sprintf("%.3f (%.3f)", cover_est,  cover_mcse)
    ) %>%
    select(all_of(grouping_vars), Method = method,
           Bias, `Rel. Bias`, `Emp. SE`, `Model SE`, Coverage)
}

plot_zip_with_estimates <- function(data, true_value) {
  data$model <- data$dgm; data$dgm <- pretty_dgm(data$dgm)
  
  data <- data %>%
    filter(!is.na(std.error) & std.error > 0) %>%
    mutate(
      z_score     = (estimate - true_value) / std.error,
      abs_z_score = abs(z_score),
      includes_true_value =
        ifelse(conf.low <= true_value & conf.high >= true_value, "Yes", "No")
    ) %>%
    group_by(method) %>%
    mutate(centile = percent_rank(abs_z_score)) %>%
    ungroup()
  
  if (nrow(data) == 0)
    return(ggplot() + labs(title = "No data.") + theme_void())
  
  perc95 <- data %>%
    group_by(method) %>%
    summarise(p95 = quantile(centile, .95), .groups = "drop")
  
  ggplot(data, aes(x = estimate, y = centile)) +
    geom_errorbarh(aes(xmin = conf.low, xmax = conf.high,
                       colour = includes_true_value), height = 0.002) +
    geom_point(colour = "white", size = 0.7) +
    geom_vline(xintercept = true_value, colour = "red", linetype = "dashed") +
    scale_color_manual(values = c(Yes = "black", No = "orange")) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    labs(x = "Estimate", y = "Fractional centile of |Z-score|",
         colour = "CI includes\ntrue value") +
    theme_minimal() +
    facet_grid(dgm ~ method) +
    geom_hline(data = perc95, aes(yintercept = p95),
               colour = "blue", linetype = "dotted") +
    ylim(0.5, 1)
}

balance_summary_tbl <- function(df_energy, df_cors) {
  ess_tbl <- df_energy %>%
    rename(eps_A = energy_A, eps_C = energy_X) %>%
    group_by(dgm, method) %>%
    summarise(
      Mean_ESS  = mean(ess), SD_ESS = sd(ess),
      P5_ESS    = quantile(ess, .05), P95_ESS = quantile(ess, .95),
      Dw_mean   = mean(D_w), Dw_sd  = sd(D_w),
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
      P95_rho  = quantile(mean_abs_cor, .95), Max_rho = max(mean_abs_cor),
      .groups = "drop"
    )
  
  ess_tbl %>%
    inner_join(cor_tbl, by = c("dgm", "method")) %>%
    select(dgm, Method = method,
           Mean_ESS, SD_ESS, P5_ESS, P95_ESS,
           Dw_mean, Dw_sd, eps_A_mean, eps_C_mean,
           Mean_rho, SD_rho, P95_rho, Max_rho) %>%
    arrange(dgm, Method)
}

round_numeric_cols <- function(df) {
  num_cols <- names(df)[vapply(df, is.numeric, logical(1))]
  for (cl in num_cols) {
    x  <- df[[cl]]; ch <- as.character(x)
    tiny  <- abs(x) < 1e-4 & x != 0
    small <- abs(x) < 1e-3 & abs(x) >= 1e-4
    ch[tiny]  <- "<0.0001"; ch[small] <- "<0.001"
    keep <- !(tiny | small)
    if (any(keep)) {
      if (grepl("ESS", cl)) {
        ch[keep] <- format(round(x[keep]), scientific = FALSE, trim = TRUE)
      } else {
        rng <- range(x[keep])
        ch[keep] <- if (max(abs(rng)) < 1e3 && max(abs(rng)) > 1e-3)
          sprintf("%.3f", round(x[keep], 3))
        else
          format(signif(x[keep], 3), scientific = TRUE)
      }
    }
    df[[cl]] <- ch
  }
  df
}

prettify_headers <- function(df) {
  nm <- names(df)
  map <- c(
    Mean_ESS  = "Mean ESS",
    SD_ESS    = "SD ESS",
    P5_ESS    = "5<sup>th</sup> perc ESS",
    P95_ESS   = "95<sup>th</sup> perc ESS",
    Dw_mean   = "Dw",
    Dw_sd     = "SD Dw",
    eps_A_mean= "&epsilon;<sub>A</sub>",
    eps_C_mean= "&epsilon;<sub>C</sub>",
    Mean_rho  = "Mean |&rho;<sub>w</sub>|",
    SD_rho    = "SD |&rho;<sub>w</sub>|",
    P95_rho   = "95<sup>th</sup> |&rho;<sub>w</sub>|",
    Max_rho   = "Max |&rho;<sub>w</sub>|"
  )
  names(df) <- ifelse(nm %in% names(map), map[nm], nm)
  df
}

# ─────────────────────────────────────────────────────────────────────────────
# 5) UI
# ─────────────────────────────────────────────────────────────────────────────
ui <- fluidPage(
  theme = bs_theme(bootswatch = "flatly", base_font = font_google("Roboto")),
  tags$style(HTML("
    div.dataTables_filter{float:left!important;text-align:left!important;}
    table.dataTable tbody td{padding:6px 10px;}
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
      checkboxGroupInput("dgm", "Exposure Distribution:", choices = NULL),
      hr(), h4("View Selection"),
      selectInput("view_type", "Result Type to View:",
                  choices = c("Model Performance" = "models",
                              "Covariate Balance" = "balance")),
      hr(), uiOutput("method_ui"),
      hr(), h4("Export"),
      uiOutput("export_ui")   # hidden in browser
    ),
    mainPanel(
      width = 9,
      tabsetPanel(id = "main_tabs",
                  tabPanel("Summary Table", DTOutput("summary_tbl")),
                  tabPanel("Zip Plot", plotOutput("zip_plot", height = "800px"))
      )
    )
  )
)

# ─────────────────────────────────────────────────────────────────────────────
# 6) Server
# ─────────────────────────────────────────────────────────────────────────────
server <- function(input, output, session) {
  
  updateSelectInput(session, "ate_exp",
                    choices = setNames(sort(unique(sim_params$ate_exp)),
                                       paste0("ln ", sort(unique(sim_params$ate_exp)))),
                    selected = 1.1
  )
  updateSelectInput(session, "mech", choices = unique(sim_params$mech))
  updateCheckboxGroupInput(session, "dgm",
                           choices  = setNames(unique(sim_params$dgm), pretty_dgm(unique(sim_params$dgm))),
                           selected = unique(sim_params$dgm)
  )
  
  output$phi_ui <- renderUI({
    req(input$mech)
    if (input$mech == "complete")
      selectInput("phi", "Proportion of Missingness:", choices = 0, selected = 0)
    else
      selectInput("phi", "Proportion of Missingness:",
                  choices = sort(unique(sim_params$phi[sim_params$mech == input$mech])))
  })
  
  output$method_ui <- renderUI({
    req(input$view_type)
    available <- if (input$view_type == "models")
      unique(sim_params$method)
    else
      setdiff(unique(sim_params$method), c("unadjusted", "adjusted"))
    checkboxGroupInput("method", "Method:", choices = available, selected = available)
  })
  
  # Show export buttons only when running locally
  output$export_ui <- renderUI({
    if (is_browser) return(NULL)
    tagList(
      downloadButton("dl_table", "Download table (Word)"),
      conditionalPanel("input.view_type == 'models'",
                       downloadButton("dl_plot", "Download plot (PDF)"))
    )
  })
  
  # Reusable filter that works for Arrow *or* data.frame
  base_filter <- function(ds) {
    ds %>%
      dplyr::filter(
        ate_exp == input$ate_exp,
        mech    == input$mech,
        phi     == input$phi,
        weight  == input$weight_type,
        dgm     %in% input$dgm,
        method  %in% input$method
      )
  }
  
  filt_energy <- reactive({
    df <- base_filter(ds_energy)
    collect_df(df) %>% mutate(dgm = pretty_dgm(dgm))
  })
  filt_cors <- reactive({
    df <- base_filter(ds_cors)
    collect_df(df) %>% mutate(dgm = pretty_dgm(dgm))
  })
  filt_models <- reactive({
    df <- base_filter(ds_models)
    collect_df(df) %>% mutate(dgm = pretty_dgm(dgm))
  })
  
  current_table_raw <- reactive({
    if (input$view_type == "models") {
      df <- filt_models(); if (!nrow(df)) return(NULL)
      sims <- rsimsum::simsum(df, "estimate", "std.error",
                              true = log(as.numeric(input$ate_exp)),
                              methodvar = "method", by = "dgm")
      format_rsimsum_summary(sims)
    } else {
      de <- filt_energy(); dc <- filt_cors()
      if (!nrow(de) || !nrow(dc)) return(NULL)
      balance_summary_tbl(de, dc)
    }
  })
  
  current_table_disp <- reactive({
    tbl <- current_table_raw(); if (is.null(tbl)) return(NULL)
    tbl <- round_numeric_cols(tbl)
    if (input$view_type == "balance") tbl <- prettify_headers(tbl)
    tbl
  })
  
  current_plot <- reactive({
    req(input$view_type == "models")
    df <- filt_models(); if (!nrow(df)) return(NULL)
    plot_zip_with_estimates(df, true_value = log(as.numeric(input$ate_exp)))
  })
  
  output$summary_tbl <- renderDT({
    if (is.null(current_table_disp()))
      return(datatable(data.frame(Message = "No data for this combination."),
                       class = "compact", rownames = FALSE))
    datatable(
      current_table_disp(),
      class = "stripe hover compact order-column row-border",
      rownames = FALSE, filter = "top",
      options = list(pageLength = 20, autoWidth = TRUE, scrollX = TRUE),
      selection = "none",
      escape = FALSE
    )
  })
  
  output$zip_plot <- renderPlot({
    g <- current_plot(); req(!is.null(g)); print(g)
  })
  observe({
    if (input$view_type != "models") hideTab("main_tabs", "Zip Plot")
    else showTab("main_tabs", "Zip Plot")
  })
  
  # Local-only downloads (not shown in browser)
  if (!is_browser && requireNamespace("officer", quietly = TRUE) &&
      requireNamespace("flextable", quietly = TRUE)) {
    output$dl_table <- downloadHandler(
      filename = function() paste0("simulation_table_", input$view_type, ".docx"),
      content = function(file) {
        tbl <- current_table_raw(); validate(need(!is.null(tbl), "No data"))
        doc <- officer::read_docx() %>%
          officer::body_add_par("Simulation summary", style = "heading 1") %>%
          officer::body_add_flextable(
            flextable::qflextable(round_numeric_cols(tbl))
          )
        print(doc, target = file)
      }
    )
    output$dl_plot <- downloadHandler(
      filename = function() "zip_plot.pdf",
      content = function(file) {
        g <- current_plot(); validate(need(!is.null(g), "No plot"))
        pdf(file, width = 8, height = 6); print(g); dev.off()
      }
    )
  }
}

# ─────────────────────────────────────────────────────────────────────────────
# 7) Run app
# ─────────────────────────────────────────────────────────────────────────────
shinyApp(ui, server)
