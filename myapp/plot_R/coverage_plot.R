# Coverage (with 95% Monte Carlo CIs) by method and exposure.
plot_coverage_by_method <- function(data, true_value, nominal = 0.95) {
  if (!nrow(data)) {
    return(ggplot2::ggplot() +
             ggplot2::labs(title = "No data.") +
             ggplot2::theme_void())
  }
  
  # Optional: pick up pretty_dgm() if defined in app.R
  .pretty <- get0("pretty_dgm", mode = "function", inherits = TRUE)
  
  # Local, NA-proof display labels for facet strips
  label_dgm_display <- function(x) {
    x <- as.character(x)
    out <- x
    out[x %in% c("nb_bin", "NegBin")]     <- "DGM: Negative Binomial"
    out[x %in% c("pois_bin", "Poisson")]  <- "DGM: Poisson"
    out[is.na(x)] <- ""
    out
  }
  
  data <- data %>%
    dplyr::mutate(
      covered = ifelse(conf.low <= true_value & conf.high >= true_value, 1, 0),
      dgm     = if (is.function(.pretty)) .pretty(dgm) else dgm,
      dgm_lab = label_dgm_display(dgm),
      method  = ifelse(method == "multinom", "multinomial", method)
    )
  
  cov_summary <- data %>%
    dplyr::group_by(method, dgm_lab) %>%
    dplyr::summarise(
      n        = dplyr::n(),
      coverage = mean(covered),
      se       = sqrt(pmax(coverage * (1 - coverage) / n, 0)),
      ci_lower = pmax(0, coverage - stats::qnorm(0.975) * se),
      ci_upper = pmin(1, coverage + stats::qnorm(0.975) * se),
      .groups  = "drop"
    )
  
  # Order methods by worst absolute deviation from nominal
  ord <- cov_summary %>%
    dplyr::group_by(method) %>%
    dplyr::summarise(worst = max(abs(coverage - nominal), na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(worst)) %>%
    dplyr::pull(method)
  cov_summary$method <- factor(cov_summary$method, levels = ord)
  
  # Data-driven y-limits
  y_min <- min(cov_summary$ci_lower, nominal, na.rm = TRUE)
  y_max <- max(cov_summary$ci_upper, nominal, na.rm = TRUE)
  pad   <- 0.02
  lims  <- c(max(0, y_min - pad), min(1, y_max + pad))
  
  # 0.05 major breaks
  brk_start <- floor(lims[1] * 20) / 20
  brk_end   <- ceiling(lims[2] * 20) / 20
  brks      <- seq(brk_start, brk_end, by = 0.05)
  
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
      legend.position  = "none",
      strip.background = ggplot2::element_rect(fill = "white", colour = NA),
      strip.text       = ggplot2::element_text(face = "bold", size = 15),
      axis.text        = ggplot2::element_text(size = 12),
      axis.title       = ggplot2::element_text(size = 13, margin = ggplot2::margin(t = 6, r = 6)),
      panel.spacing.x  = grid::unit(1.2, "lines")
    )
}
