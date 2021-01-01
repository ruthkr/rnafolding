#' Read PostScript Dot-Plot File
#'
#' @param filename The name of the PostScript file which dot-plot data is to be read from.
#'
#' @return Dot-plot data frame.
read_ps_dotplot <- function(filename) {
  # Load PostScript file
  ps_file <- readLines(filename)

  # Detect start of data
  data_start <- which(ps_file == "%start of base pair probability data")

  # Read data into data frame
  lines <- ps_file %>%
    .[data_start:length(ps_file)] %>%
    stringr::str_subset("ubox$|lbox$")

  # Detect if there is ubox data
  if (length(lines) > 0) {
    dotplot_data <- lines %>%
      textConnection() %>%
      utils::read.csv(sep = " ", header = FALSE) %>%
      `colnames<-`(c("pos_i", "pos_j", "prob", "box"))

    return(dotplot_data)
  } else {
    return(NULL)
  }
}

#' Get Position Split Intervals Data
#'
#' Get position split intervals data frame when using facets in
#' \code{plot_summary_map()} function.
#'
#' @param nt_num Number of nucleotides of the sequence.
#' @param num_facets Number of facets or splits.
#'
#' @return Split intervals data frame.
get_split_intervals <- function(nt_num, num_facets) {
  facet_width <- ceiling(nt_num / num_facets)
  intervals <- seq(0, num_facets * facet_width, facet_width)

  intervals_data <- data.frame(end = intervals) %>%
    dplyr::mutate(start = dplyr::lag(end) + 1) %>%
    dplyr::slice(-1) %>%
    tibble::rowid_to_column(var = "group") %>%
    dplyr::mutate(group = as.factor(group)) %>%
    dplyr::select(group, split_start = start, split_end = end)

  return(intervals_data)
}

#' Add Split Group to Positional Data
#'
#' Adds position split \code{group} variable to data frames when using facets
#' in \code{plot_summary_map()} function.
#'
#' @param data Data frame with \code{position} variable.
#' @param nt_num Number of nucleotides of the sequence.
#' @param num_facets Number of facets or splits.
#'
#' @return Data frame with additional grouping variable.
add_position_split_group <- function(data, nt_num, num_facets) {
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

#' Add Split Group to Region Data
#'
#' Adds position split \code{group} variable to data frames when using facets
#' in \code{plot_summary_map()} function.
#'
#' @param data Data frame describing ORFs or UTRs.
#' @param splits Split intervals data frame, result of \code{get_split_intervals()}.
#' @param type Type of region data. Either "orf" or "utr".
#'
#' @return Data frame with additional grouping variable(s).
add_region_split_group <- function(data, splits, type = c("orf", "utr")) {
  type <- match.arg(type)

  data <- data %>%
    tidyr::crossing(splits) %>%
    # Check intersection
    dplyr::mutate(
      intersect = !(split_start > end | start > split_end)
    ) %>%
    dplyr::filter(intersect) %>%
    # Fix limits
    dplyr::mutate(
      start = ifelse(start < split_start, split_start, start),
      end = ifelse(end > split_end, split_end, end)
    )

  # Select used variables
  if (type == "orf") {
    data <- data %>%
      dplyr::mutate(interval_group = as.factor(interval_group)) %>%
      dplyr::select(start, end, interval_group, group)
  } else if (type == "utr") {
    data <- data %>%
      dplyr::select(start, end, utr, group)
  }

  return(data)
}
