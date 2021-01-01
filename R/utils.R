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

#' Add Split Group to Data
#'
#' Adds position split \code{group} variable to data frames when using facets
#' in \code{plot_summary_map()} function.
#'
#' @param data Data frame with \code{position} variable.
#' @param nt_num Number of nucleotides of the sequence.
#' @param num_facets Number of facets or splits.
#'
#' @return Data frame with additional \code{group} variable.
add_split_group <- function(data, nt_num, num_facets) {
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
