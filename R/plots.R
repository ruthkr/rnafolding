#' @noRd
hide_x_axis <- function(gg) {
  gg <- gg +
    ggplot2::theme(
      axis.title.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank()
    )

  return(gg)
}

#' Individual Plots
#'
#' @param data Data frame for entropy, base-pair probabilities, ORFs, or UTRs.
#' @param prob_cutoff Probability cut-off.
#' @param split_lims Split limits data frame, result of \code{get_split_lims()}.
#' @param hide_x If TRUE, the \code{x} axis will not be shown.
#'
#' @return A \code{ggplot} object.
#'
#' @name plots
NULL
#> NULL

#' @rdname plots
#' @export
plot_entropy <- function(data, split_lims, hide_x = FALSE) {
  gg <- data %>%
    ggplot2::ggplot() +
    ggplot2::aes(x = position, y = entropy) +
    ggplot2::geom_area(fill = "lightcoral", color = "black") +
    ggplot2::scale_x_continuous(
      limits = c(split_lims[1], split_lims[2]),
      expand = c(0, 0),
      breaks = c(1, seq(0, 10000, 100))
    )

  if (hide_x) {
    gg <- hide_x_axis(gg)
  }

  return(gg)
}

#' @rdname plots
#' @export
plot_base_pairs <- function(data, prob_cutoff = 0.9, split_lims, hide_x = FALSE) {
  gg <- data %>%
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
      limits = c(split_lims[1], split_lims[2]),
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
    gg <- hide_x_axis(gg)
  }

  return(gg)
}

#' @rdname plots
#' @export
plot_orfs <- function(data, split_lims, hide_x = FALSE) {
  gg <- data %>%
    ggplot2::ggplot() +
    ggplot2::aes(
      xmin = start, ymin = 5,
      xmax = end, ymax = 10
    ) +
    ggplot2::geom_hline(yintercept = 7.5)

  if (nrow(data) > 0) {
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
      limits = c(split_lims[1], split_lims[2]),
      expand = c(0, 0),
      breaks = c(1, seq(0, 10000, 100))
    ) +
    ggplot2::scale_y_continuous(
      limits = c(4, 11)
    ) +
    ggplot2::facet_grid(interval_group ~ .) +
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
    gg <- hide_x_axis(gg)
  }

  return(gg)
}

#' @rdname plots
#' @export
plot_utrs <- function(data, split_lims, hide_x = FALSE) {
  gg <- data %>%
    ggplot2::ggplot() +
    ggplot2::geom_hline(yintercept = 7.5)

  if (nrow(data) > 0) {
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
      limits = c(split_lims[1], split_lims[2]),
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
    gg <- hide_x_axis(gg)
  }

  return(gg)
}


#' Plot Structure Map from Sliding Windows
#'
#' Plot structure map from sliding windows folding results.
#'
#' @param windowed_folds List of \code{RNAfold} results for each sliding window. Result of \code{fold()} function.
#' @param utr5_lims Vector of lower and upper limits of UTR 5'.
#' @param utr3_lims Vector of lower and upper limits of UTR 3'.
#' @param num_facets Number of facets or splits.
#' @param freq_cutoff Frequency cut-off.
#' @param prob_cutoff Probability cut-off.
#' @param plot_list List of plots to display. Options are \code{c("entropy", "bpp", "orf", "utr")}, the order of which will be respected.
#'
#' @return A \code{patchwork} object.
#' @export
plot_structure_map <- function(windowed_folds, utr5_lims = NULL, utr3_lims = NULL, num_facets = 1, freq_cutoff = 0.5, prob_cutoff = 0.9, plot_list = c("entropy", "bpp", "orf", "utr")) {
  # Sequence properties
  sequence_length <- attr(windowed_folds, "seq_length")
  seq <- attr(windowed_folds, "seq")

  # Parse plot selection
  entropy_bool <- "entropy" %in% plot_list
  bpp_bool <- "bpp" %in% plot_list
  orf_bool <- "orf" %in% plot_list
  utr_bool <- "utr" %in% plot_list
  last_plot <- plot_list[length(plot_list)]

  # List of splits
  splits_data <- get_split_lims(sequence_length, num_facets)

  # Calculate entropy
  if (entropy_bool) {
    entropy_data <- windowed_folds %>%
      aggregate_shannon_entropy(prob_cutoff = 0) %>%
      add_position_split_group(nt_num = sequence_length, num_facets)
  }

  # Calculate base-pairs
  if (bpp_bool) {
    base_pair_data <- windowed_folds %>%
      aggregate_base_pairs_prob(freq_cutoff, prob_cutoff) %>%
      # Calculate
      get_base_pairs_arcs() %>%
      add_position_split_group(nt_num = sequence_length, num_facets)
  }

  # Calculate ORFs
  if (orf_bool) {
    orf_data <- get_ORFs(seq, calc_groups = TRUE) %>%
      add_region_split_group(splits_data, type = "orf")
  }

  # Parse UTRs
  if (utr_bool) {
    utr_data <- build_UTRs_data(utr5_lims, utr3_lims) %>%
      add_region_split_group(splits_data, type = "utr")
  }

  # Create list of plots
  plot_list <- list()

  for (facet in 1:num_facets) {
    split_lims <- splits_data %>%
      dplyr::filter(group == facet) %>%
      dplyr::select(split_start, split_end) %>%
      unlist()

    if (entropy_bool) {
      plot_list[[paste0("entropy_", facet)]] <- entropy_data %>%
        dplyr::filter(group == facet) %>%
        plot_entropy(
          split_lims,
          hide_x = ("entropy" != last_plot)
        )
    }

    if (bpp_bool) {
      plot_list[[paste0("base_pair_", facet)]] <- base_pair_data %>%
        dplyr::filter(group == facet) %>%
        plot_base_pairs(
          prob_cutoff = prob_cutoff,
          split_lims,
          hide_x = ("bpp" != last_plot)
        )
    }

    if (orf_bool) {
      plot_list[[paste0("orf_", facet)]] <- orf_data %>%
        dplyr::filter(group == facet) %>%
        plot_orfs(
          split_lims,
          hide_x = ("orf" != last_plot)
        )
    }

    if (utr_bool) {
      plot_list[[paste0("utr_", facet)]] <- utr_data %>%
        dplyr::filter(group == facet) %>%
        plot_utrs(
          split_lims,
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
