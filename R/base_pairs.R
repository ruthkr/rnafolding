#' Count Window Appearances
#'
#' Count on how many sliding windows the pair \code{(i, j)} was slid over.
#'
#' @param pos_i Position of nucleotide \code{i}.
#' @param pos_j Position of nucleotide \code{j}.
#' @param windows_list List of \code{RNAfold} results for each sliding window. Result of \code{fold()} function.
#'
#' @return Integer of window appearance count.
count_window_appearances <- function(pos_i, pos_j, windows_list) {
  count_appearances <- windows_list %>%
    purrr::map(
      function(x) {
        pos_i >= x[["start_nt"]] & pos_i <= x[["end_nt"]] &
          pos_j >= x[["start_nt"]] & pos_j <= x[["end_nt"]]
      }
    ) %>%
    unlist() %>%
    sum()

  return(count_appearances)
}

#' Define Arc Segment
#'
#' @param c Vector c(x, y) defining center of the arc.
#' @param r Radius of the arc.
#' @param angles Angles between which the arc segment is defined.
#' @param length Length of points defining the arc segment.
#'
#' @return Data frame of \code{(x, y)} points defining the arc segment.
arc_segment <- function(c, r, angles = c(0, pi), length = 100) {
  seqang <- seq(angles[1], angles[2], length = length)
  x <- c[1] + r * cos(seqang)
  y <- c[2] + r * sin(seqang)

  data <- data.frame(
    arc_x = x,
    arc_y = y
  )

  return(data)
}

#' Calculate Arc Trajectory Between Two Points
#'
#' Calculates the arc trajectory between two points on the \code{x} axis.
#'
#' @param i Position of tarting point.
#' @param j Position of ending point.
#'
#' @return Arc trajectory.
get_arc_trajectory <- function(i, j) {
  pair <- c(i, j)
  center <- c(mean(pair), 0)
  radius <- abs(diff(pair)) / 2

  arc_trajectory <- dplyr::bind_cols(
    data.frame(pos_i = i, pos_j = j),
    arc_segment(center, radius)
  )

  return(arc_trajectory)
}

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
      base_pairs_data,
      by = c("pos_i", "pos_j")
    ) %>%
    dplyr::rename(
      position = arc_x,
      value = arc_y
    )

  return(arc_trajectory)
}
