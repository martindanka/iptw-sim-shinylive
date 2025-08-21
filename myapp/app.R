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
# Keep the browser install light (no ggdist here).
install_missing_pkgs(c("rsimsum", "ggplot2", "forcats"))

# 3. Base directory for app files
app_base <- if (nzchar(Sys.getenv("SHINYLIVE"))) {
  "."
} else if (dir.exists("myapp")) {
  "myapp"
} else {
  "."
}


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

data_dir <- file.path(app_base, "data")
ds_cors <- read_csv(file.path(data_dir, "cors_small.csv.gz"), show_col_types = FALSE)
ds_energy <- read_csv(file.path(data_dir, "energy_small.csv.gz"), show_col_types = FALSE)
ds_models <- read_csv(file.path(data_dir, "models_small.csv.gz"), show_col_types = FALSE)

collect_df <- identity

# Helpers & plotting ----------------------------------------------------------------------------------------------

app_theme <- bs_theme(version = 5, preset = "bootstrap")

# Idempotent mapping: accepts raw ("nb_bin","pois_bin") *or* already-pretty ("NegBin","Poisson")
pretty_dgm <- function(x) dplyr::recode(
  as.character(x),
  nb_bin = "NegBin",
  pois_bin = "Poisson",
  NegBin = "NegBin",
  Poisson = "Poisson",
  .default = as.character(x)
)

# Display labels for facet strips (pre-labelled in the data; no ggplot labeller)
label_dgm_display <- function(x) {
  x <- as.character(x)
  out <- x
  out[x %in% c("nb_bin", "NegBin")]    <- "DGM: Negative Binomial"
  out[x %in% c("pois_bin", "Poisson")] <- "DGM: Poisson"
  out[is.na(x)] <- ""  # hide any stray NA rather than printing "NA"
  out
}

sim_params <- ds_models %>%
  distinct(ate_exp, mech, phi, dgm, method, weight) %>%
  arrange(ate_exp, mech, phi)

format_rsimsum_summary <- function(simsum_obj) {
  if (is.null(simsum_obj)) return(NULL)
  
  target_stats <- c("bias", "rbias", "empse", "modelse", "cover")
  grouping_vars <- if (is.null(simsum_obj$by)) character(0) else simsum_obj$by
  
  wide <- simsum_obj$summ %>%
    filter(stat %in% target_stats) %>%
    group_by(across(all_of(c("method", grouping_vars, "stat")))) %>%
    summarise(est = mean(est), mcse = mean(mcse), .groups = "drop") %>%
    tidyr::pivot_wider(
      id_cols = c(all_of(grouping_vars), method),
      names_from = stat, values_from = c(est, mcse),
      names_glue = "{stat}_{.value}"
    )
  
  exp_cols <- as.vector(outer(target_stats, c("est", "mcse"), paste, sep = "_"))
  wide[setdiff(exp_cols, names(wide))] <- NA_real_
  
  if ("dgm" %in% grouping_vars) wide$dgm <- pretty_dgm(wide$dgm)
  
  # Build display strings; guard relative bias when undefined (e.g., true value = 0 => NaN/NA)
  to_rel_bias <- function(est, se) {
    bad <- !is.finite(est) | !is.finite(se) | is.na(est) | is.na(se)
    out <- sprintf("%.1f%% (%.1f%%)", est * 100, se * 100)
    out[bad] <- "\u2013"  # en dash only
    out
  }
  
  wide %>%
    mutate(
      `Bias (MCSE)` = sprintf("%.3f (%.3f)", bias_est, bias_mcse),
      `Rel. Bias (MCSE)` = to_rel_bias(rbias_est, rbias_mcse),
      `Emp. SE` = sprintf("%.3f", empse_est),
      `Model SE` = sprintf("%.3f", modelse_est),
      `Coverage (MCSE)` = sprintf("%.3f (%.3f)", cover_est, cover_mcse)
    ) %>%
    select(
      all_of(grouping_vars),
      Method = method,
      `Bias (MCSE)`, `Rel. Bias (MCSE)`, `Emp. SE`, `Model SE`, `Coverage (MCSE)`
    ) %>%
    rename(Exposure = dgm)
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
    return(ggplot() + labs(title = "No data.") + theme_void())
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
      title = NULL,
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
    select(
      dgm,
      Method = method,
      Mean_ESS, SD_ESS, P5_ESS, P95_ESS,
      Dw_mean, Dw_sd, eps_A_mean, eps_C_mean,
      Mean_rho, SD_rho, P95_rho, Max_rho
    ) %>%
    arrange(dgm, Method) %>%
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
      " padding:0 0.4em;",
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

# ----- Coverage plot (NA-proof facet titles, no ggplot labeller) ------------
plot_coverage_by_method <- function(data, true_value, nominal = 0.95) {
  if (!nrow(data)) {
    return(ggplot2::ggplot() +
             ggplot2::labs(title = "No data.") +
             ggplot2::theme_void())
  }
  
  .pretty <- get0("pretty_dgm", mode = "function", inherits = TRUE)
  
  data <- data %>%
    dplyr::mutate(
      covered = ifelse(conf.low <= true_value & conf.high >= true_value, 1, 0),
      dgm = if (is.function(.pretty)) .pretty(dgm) else dgm,
      dgm_lab = label_dgm_display(dgm),  # pre-labelled strip text
      method = ifelse(method == "multinom", "multinomial", method)
    )
  
  cov_summary <- data %>%
    dplyr::group_by(method, dgm_lab) %>%
    dplyr::summarise(
      n = dplyr::n(),
      coverage = mean(covered),
      se = sqrt(pmax(coverage * (1 - coverage) / n, 0)),
      ci_lower = pmax(0, coverage - stats::qnorm(0.975) * se),
      ci_upper = pmin(1, coverage + stats::qnorm(0.975) * se),
      .groups = "drop"
    )
  
  # Order methods by worst absolute deviation from nominal
  ord <- cov_summary %>%
    dplyr::group_by(method) %>%
    dplyr::summarise(worst = max(abs(coverage - nominal), na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(worst)) %>%
    dplyr::pull(method)
  cov_summary$method <- factor(cov_summary$method, levels = ord)
  
  # Data-driven y-limits; pad slightly and clip to [0,1]
  y_min <- min(cov_summary$ci_lower, nominal, na.rm = TRUE)
  y_max <- max(cov_summary$ci_upper, nominal, na.rm = TRUE)
  pad <- 0.02
  lims <- c(max(0, y_min - pad), min(1, y_max + pad))
  
  # Major breaks every 0.05; no minor breaks
  brk_start <- floor(lims[1] * 20) / 20
  brk_end <- ceiling(lims[2] * 20) / 20
  brks <- seq(brk_start, brk_end, by = 0.05)
  
  ggplot2::ggplot(cov_summary, ggplot2::aes(x = method, y = coverage)) +
    ggplot2::geom_point(size = 2, colour = "black") +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = ci_lower, ymax = ci_upper),
      width = 0, colour = "black"
    ) +
    ggplot2::geom_hline(yintercept = nominal, linetype = "dashed", linewidth = 0.6) +
    ggplot2::scale_y_continuous(
      limits = lims,
      breaks = brks,
      minor_breaks = NULL,
      expand = ggplot2::expansion(mult = c(0, 0))
    ) +
    ggplot2::labs(
      title = NULL,
      x = "Method for Constructing Weights",
      y = "Coverage with Monte Carlo 95% CI"
    ) +
    ggplot2::facet_wrap(~ dgm_lab) +
    ggplot2::theme_bw(base_size = 14) +
    ggplot2::theme(
      legend.position = "none",
      strip.background = ggplot2::element_rect(fill = "white", colour = NA),
      strip.text = ggplot2::element_text(face = "bold", size = 15),
      axis.text = ggplot2::element_text(size = 12),
      axis.title = ggplot2::element_text(size = 13, margin = ggplot2::margin(t = 6, r = 6)),
      panel.spacing.x = grid::unit(1.2, "lines")
    )
}

# ----- Bias plot (NA-proof facet titles, no ggplot labeller) ----------------
plot_bias_by_method <- function(data, true_value) {
  if (!nrow(data)) {
    return(ggplot2::ggplot() +
             ggplot2::labs(title = "No data.") +
             ggplot2::theme_void())
  }
  
  .pretty <- get0("pretty_dgm", mode = "function", inherits = TRUE)
  if (!is.function(.pretty)) .pretty <- base::identity
  
  has_ggdist <- requireNamespace("ggdist", quietly = TRUE)
  
  perf <- data %>%
    dplyr::filter(!is.na(estimate)) %>%
    dplyr::mutate(
      diff = estimate - true_value,
      method = ifelse(method == "multinom", "multinomial", method),
      dgm = .pretty(dgm),
      dgm_lab = label_dgm_display(dgm)  # pre-labelled strip text
    )
  
  bias_summary <- perf %>%
    dplyr::group_by(method, dgm_lab) %>%
    dplyr::summarise(
      n = dplyr::n(),
      bias = mean(diff),
      se_mean = stats::sd(diff) / sqrt(n),
      ci_lower = bias - stats::qnorm(0.975) * se_mean,
      ci_upper = bias + stats::qnorm(0.975) * se_mean,
      .groups = "drop"
    )
  
  # Order by |bias| (desc), then reverse so "worst" appears at TOP
  perf_ord <- perf %>%
    dplyr::left_join(bias_summary, by = c("method", "dgm_lab")) %>%
    dplyr::mutate(
      method = forcats::fct_reorder(method, abs(bias), .desc = TRUE),
      method = if ("adjusted" %in% levels(method))
        forcats::fct_relevel(method, "adjusted", after = 6) else method
    )
  lev_rev <- rev(levels(perf_ord$method))
  perf_ord$method <- factor(perf_ord$method, levels = lev_rev)
  bias_summary$method <- factor(bias_summary$method, levels = lev_rev)
  
  p <- ggplot2::ggplot(perf_ord, ggplot2::aes(x = method, y = diff, fill = method))
  
  if (has_ggdist) {
    p <- p + ggdist::stat_halfeye(
      adjust = 0.5, width = 0.6, justification = -0.2, .width = 0, point_colour = NA
    )
  } else {
    p <- p + ggplot2::geom_violin(width = 0.6, trim = TRUE, alpha = 0.9, colour = NA)
  }
  
  p +
    ggplot2::geom_point(
      data = bias_summary,
      ggplot2::aes(x = method, y = bias),
      colour = "black", size = 2, inherit.aes = FALSE
    ) +
    ggplot2::geom_segment(
      data = bias_summary,
      ggplot2::aes(x = as.numeric(method), xend = as.numeric(method), y = 0, yend = bias),
      colour = "black", inherit.aes = FALSE
    ) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.6) +
    ggplot2::scale_y_continuous(
      limits = c(-0.5, 0.5),
      breaks = seq(-0.5, 0.5, 0.1),
      minor_breaks = NULL
    ) +
    ggplot2::labs(
      title = NULL,
      x = "Method for Constructing Weights",
      y = expression("Bias on the Log Risk Ratio Scale and Distribution of" ~
                       hat(italic(theta)[i]) - italic(theta))
    ) +
    ggplot2::facet_wrap(~ dgm_lab) +
    ggplot2::coord_flip() +
    ggplot2::theme_bw(base_size = 14) +
    ggplot2::theme(
      legend.position = "none",
      strip.background = ggplot2::element_rect(fill = "white", colour = NA),
      strip.text = ggplot2::element_text(face = "bold", size = 15),
      axis.text = ggplot2::element_text(size = 12),
      axis.title = ggplot2::element_text(size = 13, margin = ggplot2::margin(t = 6, r = 6)),
      panel.grid = ggplot2::element_blank(),
      panel.spacing.x = grid::unit(1.2, "lines")
    )
}

# --------------------------------- UI --------------------------------------

ui <- fluidPage(
  ## Roboto at run-time (no build-time effect) -------------------------------
  tags$head(
    tags$link(
      rel = "stylesheet",
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
        var img = findVisiblePlotImage(['bias_plot','coverage_plot','zip_plot']) ||
                  (document.getElementById('bias_plot')||{}).getElementsByTagName?.('img')[0] ||
                  (document.getElementById('coverage_plot')||{}).getElementsByTagName?.('img')[0] ||
                  (document.getElementById('zip_plot')||{}).getElementsByTagName?.('img')[0];

        if (!img) { alert('Plot image not found'); return; }
        var doc = new jsPDF();
        doc.addImage(img.src, 'PNG', 10, 10, 190, 0);
        doc.save(msg.filename);
      });
    "))
  ),
  
  theme = app_theme,
  
  tags$style(HTML("
  div.dataTables_filter{ float:left!important; text-align:left!important; }
  table.dataTable tbody td{ padding:6px 10px; }
  table.dataTable thead th{
    text-align:center;
    vertical-align:bottom;
    padding:6px 16px 4px 16px;
    white-space:nowrap;
  }
")),
  
  titlePanel("Simulation Results for IPW of Count Exposures"),
  
  sidebarLayout(
    sidebarPanel(
      width = 3,
      
      h4("Simulation Parameters"),
      selectInput("ate_exp", "Effect Size (ATE):", choices = NULL),
      uiOutput("mech_ui"),
      uiOutput("phi_ui"),
      selectInput(
        "weight_type", "Winsorise at 99th perc:",
        choices = c("No" = "raw", "Yes" = "win"), selected = "raw"
      ),
      checkboxGroupInput("dgm", "Exposure Distribution:", choices = NULL),
      
      hr(),
      h4("View Selection"),
      selectInput(
        "view_type", "Result Type to View:",
        choices = c("Model Performance" = "models",
                    "Covariate Balance" = "balance")
      ),
      
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
        tabPanel("Bias Plot", plotOutput("bias_plot", height = "800px")),
        tabPanel("Coverage Plot", plotOutput("coverage_plot", height = "800px")),
        tabPanel("Zip Plot", plotOutput("zip_plot", height = "800px"))
      )
    )
  )
)


# -------------------------------- Server -----------------------------------

server <- function(input, output, session) {
  ## Populate selectors -------------------------------------------------------
  updateSelectInput(
    session, "ate_exp",
    choices = setNames(sort(unique(sim_params$ate_exp)),
                       paste0("ln ", sort(unique(sim_params$ate_exp)))),
    selected = 1.1
  )
  # Missingness mechanism appears conditionally via output$mech_ui (see below)
  updateCheckboxGroupInput(
    session, "dgm",
    choices = setNames(unique(sim_params$dgm), pretty_dgm(unique(sim_params$dgm))),
    selected = "nb_bin"
  )
  
  # Conditional UI for missingness mechanism (only when ATE == 1.1)
  output$mech_ui <- renderUI({
    req(input$ate_exp)
    if (isTRUE(as.numeric(input$ate_exp) == 1.1)) {
      selectInput("mech", "Missingness Mechanism:", choices = unique(sim_params$mech))
    } else {
      NULL
    }
  })
  
  output$phi_ui <- renderUI({
    req(input$ate_exp)
    if (!isTRUE(as.numeric(input$ate_exp) == 1.1)) return(NULL)
    req(input$mech)
    if (input$mech == "complete") {
      selectInput("phi", "Proportion of Missingness:", choices = 0, selected = 0)
    } else {
      selectInput(
        "phi", "Proportion of Missingness:",
        choices = sort(unique(sim_params$phi[sim_params$mech == input$mech]))
      )
    }
  })
  
  output$method_ui <- renderUI({
    req(input$view_type)
    available <- if (input$view_type == "models") {
      unique(sim_params$method)
    } else {
      setdiff(unique(sim_params$method), c("unadjusted", "adjusted"))
    }
    available_sorted <- available[order(tolower(available))]
    checkboxGroupInput("method", "Method:", choices = available_sorted, selected = available_sorted)
  })
  
  ## Reactive helpers for mech/phi so filters behave when controls are hidden --
  mech_value <- reactive({
    if (isTRUE(as.numeric(input$ate_exp) == 1.1)) {
      req(input$mech)
    } else {
      "complete"
    }
  })
  phi_value <- reactive({
    if (!isTRUE(as.numeric(input$ate_exp) == 1.1)) {
      0
    } else {
      req(input$mech)
      if (input$mech == "complete") 0 else req(input$phi)
    }
  })
  
  ## Reactive filters ---------------------------------------------------------
  base_filter <- function(ds) {
    ds %>% filter(
      ate_exp == input$ate_exp,
      mech == mech_value(),
      phi == phi_value(),
      weight == input$weight_type,
      dgm %in% input$dgm,
      method %in% input$method
    )
  }
  filt_energy <- reactive(collect_df(base_filter(ds_energy)) %>% mutate(dgm = pretty_dgm(dgm)))
  filt_cors <- reactive(collect_df(base_filter(ds_cors)) %>% mutate(dgm = pretty_dgm(dgm)))
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
    # Add dagger to the header when true value = 0 (ATE = 1.0)
    if (input$view_type == "models" && isTRUE(as.numeric(input$ate_exp) == 1.0)) {
      names(tbl)[names(tbl) == "Rel. Bias (MCSE)"] <- "Rel. Bias (MCSE)<sup>\u2020</sup>"
    }
    if (input$view_type == "balance") tbl <- prettify_headers(tbl)
    tbl
  })
  
  output$summary_tbl <- renderDT({
    tbl <- current_table_disp()
    if (is.null(tbl)) {
      datatable(
        data.frame(Message = "No data for this combination."),
        class = "compact",
        rownames = FALSE
      )
    } else {
      
      ## ---------- FOOTNOTE TEXT & CAPTION (updated) ----------------------- ##
      footnote <- if (input$view_type == "models") {
        base_ft <- paste(
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
        if (isTRUE(as.numeric(input$ate_exp) == 1.0)) {
          paste0(
            base_ft,
            "<br /><span style='font-style:italic'>\u2020 Relative bias is undefined when the true value is zero (RR = 1.0), so cells show a dash.</span>"
          )
        } else {
          base_ft
        }
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
      }
      ## -------------------------------------------------------------------- ##
      
      # Build DataTables options with centred cells except Exposure & Method;
      # and left-align the headers of those two columns ONLY on the models view.
      ncols <- ncol(tbl)
      all_idx0 <- 0:(ncols - 1)                     # 0-based for DataTables
      left_keep <- which(colnames(tbl) %in% c("Exposure", "Method"))
      centre_targets <- setdiff(all_idx0, left_keep - 1)
      opts <- list(
        scrollX = TRUE,
        pageLength = 20,
        autoWidth = TRUE,
        fixedColumns = list(leftColumns = if (input$view_type == "balance") 2 else 1)
      )
      if (length(centre_targets)) {
        opts$columnDefs <- list(list(className = "dt-center", targets = centre_targets))
      }
      if (input$view_type == "models") {
        # Left-align the header text for Exposure and Method
        opts$headerCallback <- DT::JS("
          function(thead, data, start, end, display) {
            var $ths = $(thead).find('th');
            if ($ths.length >= 2) {
              $ths.eq(0).css('text-align', 'left');
              $ths.eq(1).css('text-align', 'left');
            }
          }
        ")
      }
      
      datatable(
        tbl,
        class = "stripe hover compact order-column row-border",
        rownames = FALSE,
        filter = "top",
        extensions = c("FixedColumns"),
        options = opts,
        selection = "none",
        escape = FALSE,
        caption = tags$caption(
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
  
  ## Bias plot ----------------------------------------------------------------
  current_bias_plot <- reactive({
    req(input$view_type == "models")
    df <- collect_df(base_filter(ds_models))
    if (!nrow(df)) return(NULL)
    plot_bias_by_method(df, true_value = log(as.numeric(input$ate_exp)))
  })
  
  output$bias_plot <- renderPlot({
    g <- current_bias_plot(); req(!is.null(g))
    print(g)
  })
  
  observe({
    if (input$view_type != "models") {
      hideTab("main_tabs", "Bias Plot")
      hideTab("main_tabs", "Coverage Plot")
      hideTab("main_tabs", "Zip Plot")
    } else {
      showTab("main_tabs", "Bias Plot")
      showTab("main_tabs", "Coverage Plot")
      showTab("main_tabs", "Zip Plot")
    }
  })
  
  ## CLIENT-SIDE DOWNLOADS ----------------------------------------------------
  # 1) Table as CSV
  observeEvent(input$btn_dl_csv, {
    tbl <- current_table_disp()
    validate(need(!is.null(tbl), "No data to download"))
    
    # Clean header names for CSV for all views (strip HTML + entities)
    clean_headers <- function(x) {
      x <- gsub("<[^>]+>", "", x)             # strip tags (incl. <sup>†</sup>)
      x <- gsub("&epsilon;", "epsilon", x, fixed = TRUE)
      x <- gsub("&rho;", "rho", x, fixed = TRUE)
      x <- gsub("\u00A0", " ", x, fixed = TRUE)  # non-breaking space
      x <- gsub("epsilon([AC])", "epsilon_\\L\\1", x, perl = TRUE)
      x <- gsub("rho([w])", "rho_\\1", x, perl = TRUE)
      x <- gsub("([0-9])<sup>th</sup>", "\\1th", x)  # just in case any linger
      x
    }
    coln <- clean_headers(colnames(tbl))
    
    # Ensure no HTML entity for en dash in data (defensive; we use Unicode)
    tbl_csv <- tbl
    tbl_csv[] <- lapply(tbl_csv, function(col) {
      if (is.character(col)) gsub("&ndash;", "\u2013", col, fixed = TRUE) else col
    })
    
    csv_lines <- c(
      paste(coln, collapse = ","),
      apply(tbl_csv, 1, function(r) paste(r, collapse = ","))
    )
    csv_text <- paste(csv_lines, collapse = "\r\n")
    session$sendCustomMessage(
      "download-csv",
      list(
        data = csv_text,
        filename = paste0("summary-", Sys.Date(), ".csv")
      )
    )
  })
  
  # 2) Plot as PDF
  observeEvent(input$btn_dl_pdf, {
    session$sendCustomMessage(
      "download-pdf",
      list(filename = paste0("plot-", Sys.Date(), ".pdf"))
    )
  })
}

# Run app ---------------------------------------------------------------------------------------------------------

shinyApp(ui, server)
