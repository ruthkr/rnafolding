#' Aggregate Base-Pair Probabilities
#'
#' @param windowed_folds List of \code{RNAfold} results for each sliding window. Result of \code{fold()} function.
#' @param freq_cutoff Frequency cut-off.
#' @param prob_cutoff Probability cut-off.
#'
#' @return Data frame specifying folding probability of the base-pair \code{(pos_i, pos_j)}.
#' @export
aggregate_base_pairs_prob <- function(windowed_folds, freq_cutoff = 0.5, prob_cutoff = 0.9) {
  # List of windows
  windows_list <- windowed_folds %>%
    purrr::map(
      function(x) {
        c(
          start_nt = attr(x, "start_n"),
          end_nt = attr(x, "end_n")
        )
      }
    )

  # Calculate positional frequencies
  base_pairs_data <- windowed_folds %>%
    # Bind all dot-plot data
    purrr::map(
      function(x) {
        return(x$dotplot_data)
      }
    ) %>%
    purrr::reduce(dplyr::bind_rows) %>%
    # Filter probabilities
    dplyr::filter(prob >= prob_cutoff) %>%
    # Calculate frequency
    dplyr::group_by(pos_i, pos_j) %>%
    dplyr::summarise(
      pair_count = dplyr::n(),
      prob = mean(prob),
      .groups = "drop"
    ) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      count_appearances = count_window_appearances(pos_i, pos_j, windows_list),
      frequency = pair_count / count_appearances
    ) %>%
    dplyr::filter(frequency >= freq_cutoff) %>%
    dplyr::select(pos_i, pos_j, prob) %>%
    # Clean base-pairs
    clean_base_pairs() %>%
    dplyr::select(-rowid)

  return(base_pairs_data)
}

#' Clean Base-Pairs Data
#'
#' Clean base-pairs data to ensure no nucleotide is paired with more than one nucleotide.
#'
#' @param base_pairs_data Data frame specifying folding probability of the base-pair \code{(pos_i, pos_j)}.
#'
#' @return Data frame with unique base-pair foldings.
clean_base_pairs <- function(base_pairs_data) {
  # Raw data
  data <- base_pairs_data %>%
    dplyr::select(pos_i, pos_j, prob) %>%
    dplyr::distinct() %>%
    tibble::rowid_to_column()

  # Duplicate data with swapped positions
  data <- dplyr::bind_rows(
    data,
    data %>%
      dplyr::select(rowid, pos_i = pos_j, pos_j = pos_i, prob)
  ) %>%
    dplyr::ungroup()

  # Select highest probability
  data_clean <- data

  num_nt <- max(data_clean$pos_i, data_clean$pos_j)

  for (nt in 1:num_nt) {
    data_small <- data %>%
      dplyr::filter(pos_i == nt)

    if (nrow(data_small) > 0) {
      good_row <- data_small %>%
        dplyr::arrange(dplyr::desc(prob)) %>%
        dplyr::slice(1)

      bad_ids <- data_small %>%
        dplyr::filter(rowid != good_row$rowid) %>%
        dplyr::pull(rowid)

      if (length(bad_ids) > 0) {
        data_clean <- data_clean %>%
          dplyr::filter(!(rowid %in% bad_ids))
      }
    }
  }

  # Final cleanup
  data_clean <- data_clean %>%
    # Delete duplicates
    dplyr::group_by(rowid) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup()

  return(data_clean)
}

#' Calculate Base-Pairs Arc Trajectories
#'
#' Calculate base-pairs segments to build base-pairing probabilities plot.
#'
#' @param base_pairs_data Data frame specifying folding probability of the base-pair \code{(pos_i, pos_j)}.
#'
#' @return Data frame arc trajectories for each base-pair.
get_base_pairs_arcs <- function(base_pairs_data) {
  arc_trajectory <- base_pairs_data %>%
    t() %>%
    as.data.frame() %>%
    purrr::map(
      function(x) {
        get_arc_trajectory(x[1], x[2])
      }
    ) %>%
    purrr::reduce(dplyr::bind_rows)

  # Include probability
  arc_trajectory <- arc_trajectory %>%
    dplyr::left_join(
      positional_freq_data,
      by = c("pos_i", "pos_j")
    ) %>%
    dplyr::rename(
      position = arc_x,
      value = arc_y
    )

  return(arc_trajectory)
}
