#' Calculate and Aggregate Shannon's Entropy
#'
#' @param dotplot_data Dot-plot data frame.
#' @param windowed_folds List of \code{RNAfold} results for each sliding window. Result of \code{fold()} function.
#' @param prob_cutoff Probability cut-off.
#'
#' @name entropy
NULL
#> NULL

#' @rdname entropy
#' @export
calc_shannon_entropy <- function(dotplot_data, prob_cutoff = 0) {
  data <- dotplot_data %>%
    dplyr::filter(prob >= prob_cutoff) %>%
    dplyr::mutate(
      prob = prob * prob,
      entropy = ifelse(prob > 0, prob * log10(prob), 0)
    )

  entropy_data <- dplyr::bind_rows(
    data %>%
      dplyr::group_by(position = pos_i) %>%
      dplyr::summarise(
        entropy = sum(entropy),
        .groups = "drop"
      ),
    data %>%
      dplyr::group_by(position = pos_j) %>%
      dplyr::summarise(
        entropy = sum(entropy),
        .groups = "drop"
      )
  ) %>%
    # Calculate entropy
    dplyr::group_by(position) %>%
    dplyr::summarise(
      entropy = -sum(entropy),
      .groups = "drop"
    ) %>%
    dplyr::ungroup()

  return(entropy_data)
}

#' @rdname entropy
#' @export
aggregate_shannon_entropy <- function(windowed_folds, prob_cutoff = 0) {
  entropy <- windowed_folds %>%
    purrr::map(
      function(x) {
        entropy <- calc_shannon_entropy(x$dotplot_data, prob_cutoff)
        return(entropy)
      }
    ) %>%
    purrr::reduce(dplyr::bind_rows) %>%
    dplyr::group_by(position) %>%
    dplyr::summarise(
      entropy = stats::median(entropy),
      .groups = "drop"
    )

  return(entropy)
}
