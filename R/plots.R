split_plot <- function(data, nt_num, num_facets) {
  facet_width <- ceiling(nt_num / num_facets)

  data <- data %>%
    dplyr::mutate(
      group = cut(
        x = position,
        breaks = seq(0, num_facets * facet_width, facet_width)
      ),
      group = factor(group, labels = 1:num_facets)
    )

  return(data)
}

plot_entropy <- function(entropy_data, split_limits, hide_x = FALSE) {
  gg <- entropy_data %>%
    ggplot2::ggplot() +
    ggplot2::aes(x = position, y = entropy) +
    ggplot2::geom_area(fill = "lightcoral", color = "black") +
    ggplot2::scale_x_continuous(
      limits = c(split_limits[1], split_limits[2]),
      expand = c(0, 0),
      breaks = c(1, seq(0, 10000, 100))
    )

  if (hide_x) {
    gg <- gg +
      ggplot2::theme(
        # x-Axis
        axis.title.x = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank()
      )
  }

  return(gg)
}

plot_base_pairing <- function(base_pairing_data, prob_cutoff = 0.9, split_limits, hide_x = FALSE) {
  gg <- base_pairing_data %>%
    ggplot2::ggplot() +
    ggplot2::aes(
      x = position,
      y = value,
      color = prob,
      group = interaction(pos_i, pos_j)
    ) +
    ggplot2::geom_line(size = 0.5) +
    viridis::scale_color_viridis(limits = c(prob_cutoff, 1), direction = -1) +
    ggplot2::scale_x_continuous(
      limits = c(split_limits[1], split_limits[2]),
      expand = c(0, 0),
      breaks = c(1, seq(0, 10000, 100))
    ) +
    ggplot2::labs(y = "bpp") +
    ggplot2::theme(
      # y-Axis
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor.y = ggplot2::element_blank(),
    )

  if (hide_x) {
    gg <- gg +
      ggplot2::theme(
        # x-Axis
        axis.title.x = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank()
      )
  }

  return(gg)
}

plot_orfs <- function(orf_data, split_limits, hide_x = FALSE) {
  gg <- orf_data %>%
    ggplot2::ggplot() +
    ggplot2::aes(
      xmin = start, ymin = 5,
      xmax = end, ymax = 10
    ) +
    ggplot2::geom_hline(yintercept = 7.5)

  if (nrow(orf_data) > 0) {
    gg <- gg +
      ggplot2::aes(
        xmin = start, ymin = 5,
        xmax = end, ymax = 10,
        fill = utr
      ) +
      ggplot2::geom_rect(fill = "cornflowerblue", alpha = 1, color = "black")
  }

  gg <- gg +
    ggplot2::scale_x_continuous(
      limits = c(split_limits[1], split_limits[2]),
      expand = c(0, 0),
      breaks = c(1, seq(0, 10000, 100))
    ) +
    ggplot2::scale_y_continuous(
      limits = c(4, 11)
    ) +
    ggplot2::facet_grid(int_group ~ .) +
    ggplot2::labs(y = "ORFs") +
    ggplot2::theme(
      # Facet stripts
      strip.text = ggplot2::element_blank(),
      strip.text.y = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank(),
      panel.spacing = ggplot2::unit(0, "lines"),
      # y-Axis
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor.y = ggplot2::element_blank()
    )


  if (hide_x) {
    gg <- gg +
      ggplot2::theme(
        axis.title.x = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank()
      )
  }

  return(gg)
}

plot_utrs <- function(utr_data, split_limits, hide_x = FALSE) {
  gg <- utr_data %>%
    ggplot2::ggplot() +
    ggplot2::geom_hline(yintercept = 7.5)

  if (nrow(utr_data) > 0) {
    gg <- gg +
      ggplot2::aes(
        xmin = start, ymin = 5,
        xmax = end, ymax = 10,
        fill = utr
      ) +
      ggplot2::geom_rect(alpha = 1, color = "black") +
      ggplot2::scale_fill_manual(values = c("3'UTR" = "plum3", "5'UTR" = "mediumseagreen"))
  }

  gg <- gg +
    ggplot2::scale_x_continuous(
      limits = c(split_limits[1], split_limits[2]),
      expand = c(0, 0),
      breaks = c(1, seq(0, 10000, 100))
    ) +
    ggplot2::scale_y_continuous(
      limits = c(4, 11)
    ) +
    ggplot2::labs(y = "UTRs") +
    ggplot2::theme(
      # Facet stripts
      strip.text = ggplot2::element_blank(),
      strip.text.y = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank(),
      panel.spacing = ggplot2::unit(0, "lines"),
      # y-Axis
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor.y = ggplot2::element_blank(),
      legend.title = ggplot2::element_blank()
    )

  if (hide_x) {
    gg <- gg +
      ggplot2::theme(
        axis.title.x = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank()
      )
  }

  return(gg)
}

plot_structure_map <- function(window_results, utr5_lims = NULL, utr3_lims = NULL, num_facets = 3, freq_cutoff = 0.5, prob_cutoff = 0.9, plot_list = c("entropy", "bpp", "orf", "utr")) {
  # Sequence properties
  sequence_length <- attr(window_results, "seq_length")

  # Parse plot selection
  entropy_bool <- "entropy" %in% plot_list
  bpp_bool <- "bpp" %in% plot_list
  orf_bool <- "orf" %in% plot_list
  utr_bool <- "utr" %in% plot_list
  last_plot <- plot_list[length(plot_list)]

  # List of splits
  splits <- get_split_intervals(sequence_length, num_facets)

  # Calculate Entropy
  if (entropy_bool) {
    entropy_data <- window_results %>%
      aggregate_shannon_entropy(prob_cutoff = 0) %>%
      split_plot(nt_num = sequence_length, num_facets)
  }

  # Calculate Folding
  if (bpp_bool) {
    folding <- window_results %>%
      aggregate_positional_freq(freq_cutoff, prob_cutoff) %>%
      clean_folding_pairs() %>%
      dplyr::select(-rowid)

    base_pair_data <- calculate_folding_pairs(folding) %>%
      split_plot(nt_num = sequence_length, num_facets)
  }

  # Calculate ORFs
  if (orf_bool) {
    orf_data <- get_ORFs(window_results) %>%
      split_ORFs(splits)
  }

  # Parse UTRs
  if (utr_bool) {
    utr_data <- build_UTRs_data(utr5_lims, utr3_lims) %>%
      split_UTRs(splits)
  }

  # Create list of plots
  plot_list <- list()

  for (facet in 1:num_facets) {
    split_limits <- splits %>%
      dplyr::filter(group == facet) %>%
      dplyr::select(split_start, split_end) %>%
      unlist()

    if (entropy_bool) {
      plot_list[[paste0("entropy_", facet)]] <- entropy_data %>%
        dplyr::filter(group == facet) %>%
        plot_entropy(
          split_limits,
          hide_x = ("entropy" != last_plot)
        )
    }

    if (bpp_bool) {
      plot_list[[paste0("base_pair_", facet)]] <- base_pair_data %>%
        dplyr::filter(group == facet) %>%
        plot_base_pairing(
          prob_cutoff = prob_cutoff,
          split_limits,
          hide_x = ("bpp" != last_plot)
        )
    }

    if (orf_bool) {
      plot_list[[paste0("orf_", facet)]] <- orf_data %>%
        dplyr::filter(group == facet) %>%
        plot_orfs(
          split_limits,
          hide_x = ("orf" != last_plot)
        )
    }

    if (utr_bool) {
      plot_list[[paste0("utr_", facet)]] <- utr_data %>%
        dplyr::filter(group == facet) %>%
        plot_utrs(
          split_limits,
          hide_x = ("utr" != last_plot)
        )
    }
  }

  # Assemble patchwork
  patch <- plot_list %>%
    patchwork::wrap_plots() +
    patchwork::plot_layout(ncol = 1, guides = "collect")

  return(patch)
}
