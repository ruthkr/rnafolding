#' Get Open Reading Frame Data
#'
#' Find open reading frames (ORFs) of the sequence. Optionally, the ORFs can
#' be separated into non-overlapping groups, allowing to easily plot them as
#' genome map plots.
#'
#' @param seq Sequence input.
#' @param calc_groups If TRUE, the ORFs will be separated into non-overlapping groups.
#'
#' @return Data frame specifying \code{start}, \code{end}, and optionally \code{interval_group} for each ORF.
#' @export
get_ORFs <- function(seq, calc_groups = TRUE) {
  # Calculate ORFs
  orf_data <- ORFik::findORFs(seq) %>%
    as.data.frame() %>%
    dplyr::select(start, end)

  if (calc_groups) {
    # ORFs as Invervals object
    orf_intervals <- orf_data %>%
      as.matrix() %>%
      intervals::Intervals()

    # Calculate non-overlapping interval groups
    interval_groups <- intervals::detach_overlaps(orf_intervals)

    orf_data <- data.frame(
      orf_data %>%
        dplyr::select(start, end),
      interval_group = interval_groups
    )
  }

  return(orf_data)
}

#' Build UTRs Data Frame
#'
#' @param utr5_lims Vector of lower and upper limits of UTR 5'.
#' @param utr3_lims Vector of lower and upper limits of UTR 3'
#'
#' @return Data frame specifying \code{start} and \code{end} of each UTR.
build_UTRs_data <- function(utr5_lims = NULL, utr3_lims = NULL) {
  # Initialise data frame
  utr_data <- data.frame(
    utr = character(),
    start = integer(),
    end = integer()
  )

  # Bind UTR5 row
  if (length(utr5_lims) == 2) {
    utr_data <- utr_data %>%
      dplyr::bind_rows(
        list(
          utr = "5'UTR",
          start = as.integer(utr5_lims[[1]]),
          end = as.integer(utr5_lims[[2]])
        )
      )
  }

  # Bind UTR3 row
  if (length(utr3_lims) == 2) {
    utr_data <- utr_data %>%
      dplyr::bind_rows(
        list(
          utr = "3'UTR",
          start = as.integer(utr3_lims[[1]]),
          end = as.integer(utr3_lims[[2]])
        )
      )
  }

  return(utr_data)
}


