
# Coverage (with 95% Monte Carlo CIs) by method and exposure.
plot_coverage_by_method <- function(data, true_value, nominal = 0.95) {
  if (!nrow(data)) {
    return(ggplot2::ggplot() + ggplot2::labs(title = "No data.") + ggplot2::theme_void())
  }
  
  # Safely look up pretty_dgm() if it exists in the app
  .pretty <- get0("pretty_dgm", mode = "function", inherits = TRUE)
  
  data <- data %>%
    dplyr::mutate(
      covered = ifelse(conf.low <= true_value & conf.high >= true_value, 1, 0),
      dgm = if (is.function(.pretty)) .pretty(dgm) else dgm
    )
  
  cov_summary <- data %>%
    dplyr::group_by(method, dgm) %>%
    dplyr::summarise(
      n = dplyr::n(),
      coverage = mean(covered),
      se = sqrt(pmax(coverage * (1 - coverage) / n, 0)),
      ci_lower = pmax(0, coverage - qnorm(0.975) * se),
      ci_upper = pmin(1, coverage + qnorm(0.975) * se),
      .groups = "drop"
    )
  
  ggplot2::ggplot(cov_summary, ggplot2::aes(x = method, y = coverage)) +
    ggplot2::geom_point() +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = ci_lower, ymax = ci_upper), width = 0) +
    ggplot2::geom_hline(yintercept = nominal, linetype = "dashed", colour = "red") +
    ggplot2::facet_wrap(~ dgm) +
    ggplot2::labs(
      x = "Method",
      y = "Coverage",
      title = "Coverage by Method and Exposure (95% Monte Carlo CI)"
    ) +
    ggplot2::theme_bw()
}
