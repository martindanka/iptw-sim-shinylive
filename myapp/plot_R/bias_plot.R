# Half-eye ("rainbow") bias plot with 95% Monte Carlo CIs by method and exposure.

plot_bias_by_method <- function(data, true_value) {
  if (!nrow(data)) {
    return(ggplot2::ggplot() +
             ggplot2::labs(title = "No data.") +
             ggplot2::theme_void())
  }
  
  .pretty <- get0("pretty_dgm", mode = "function", inherits = TRUE)
  if (!is.function(.pretty)) .pretty <- base::identity
  
  has_ggdist <- requireNamespace("ggdist", quietly = TRUE)
  
  perf <- data |>
    dplyr::filter(!is.na(estimate)) |>
    dplyr::mutate(
      diff   = estimate - true_value,
      method = ifelse(method == "multinom", "multinomial", method),
      dgm    = .pretty(dgm)
    )
  
  bias_summary <- perf |>
    dplyr::group_by(method, dgm) |>
    dplyr::summarise(
      n        = dplyr::n(),
      bias     = mean(diff),
      se_mean  = stats::sd(diff) / sqrt(n),
      ci_lower = bias - stats::qnorm(0.975) * se_mean,
      ci_upper = bias + stats::qnorm(0.975) * se_mean,
      .groups  = "drop"
    )
  
  perf_ord <- perf |>
    dplyr::left_join(bias_summary, by = c("method", "dgm")) |>
    dplyr::mutate(
      method = forcats::fct_reorder(method, abs(bias), .desc = TRUE),
      method = if ("adjusted" %in% levels(method))
        forcats::fct_relevel(method, "adjusted", after = 6) else method
    )
  
  # Show "worst" at the TOP
  lev_rev <- rev(levels(perf_ord$method))
  perf_ord$method <- factor(perf_ord$method, levels = lev_rev)
  bias_summary <- bias_summary |>
    dplyr::mutate(method = factor(method, levels = lev_rev))
  
  lab_dgm <- ggplot2::as_labeller(c(
    "NegBin"  = "DGM: Negative Binomial",
    "Poisson" = "DGM: Poisson"
  ))
  
  p <- ggplot2::ggplot(perf_ord, ggplot2::aes(x = method, y = diff, fill = method))
  
  if (has_ggdist) {
    # Preferred look
    p <- p +
      ggdist::stat_halfeye(
        adjust = 0.5, width = 0.6, justification = -0.2, .width = 0, point_colour = NA
      )
  } else {
    # Lightweight fallback (no ggdist): narrow violin
    p <- p +
      ggplot2::geom_violin(width = 0.6, trim = TRUE, alpha = 0.9, colour = NA)
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
    ggplot2::scale_y_continuous(limits = c(-0.5, 0.5),
                                breaks = seq(-0.5, 0.5, 0.1),
                                minor_breaks = NULL) +
    ggplot2::labs(
      x = "Method for Constructing Weights",
      y = expression("Bias on the Log Risk Ratio Scale and Distribution of" ~
                       hat(italic(theta)[i]) - italic(theta))
    ) +
    ggplot2::facet_wrap(~ dgm, labeller = lab_dgm) +
    ggplot2::coord_flip() +
    ggplot2::theme_bw(base_size = 14) +
    ggplot2::theme(
      legend.position  = "none",
      strip.background = ggplot2::element_rect(fill = "white", colour = NA),
      strip.text       = ggplot2::element_text(face = "bold", size = 15),
      axis.text        = ggplot2::element_text(size = 12),
      axis.title       = ggplot2::element_text(size = 13, margin = ggplot2::margin(t = 6, r = 6)),
      panel.grid       = ggplot2::element_blank(),
      panel.spacing.x  = grid::unit(1.2, "lines")
    )
}
