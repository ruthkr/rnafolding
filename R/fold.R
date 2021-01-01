#' Fold with Sliding Windows
#'
#' @param filename
#' @param winsize
#' @param span
#' @param stepsize
#' @param increased_sample
#' @param rnafold_params
#' @param verbose
#'
#' @return
#' @export
#'
#' @examples
fold <- function(filename, winsize, span = NULL, stepsize, increased_sample = TRUE, rnafold_params = "-p", verbose = FALSE) {
  # Read whole sequence
  fasta <- seqinr::read.fasta(filename)
  seq <- fasta[[1]] %>%
    toupper()
  seq_name <- attr(fasta[[1]], "Annot") %>%
    stringr::str_remove_all(">") %>%
    stringr::str_trim()

  # Calculate number of nucleotides and windows
  num_nt <- length(seq)
  windows <- ceiling(num_nt / stepsize)

  # Calculate first window index
  if (increased_sample) {
    initial_window <- -floor(winsize / stepsize) + 1
    if (winsize %% stepsize == 0) {
      initial_window <- initial_window + 1
    }
  } else {
    initial_window <- 1
  }

  # Prepare temporary file paths
  temp_dir <- tempdir()
  temp_fasta_path <- paste0(temp_dir, "/", "temp.fasta")
  temp_seq_name <- "tempseq"
  temp_dp_path <- paste0(temp_dir, "/", temp_seq_name, "_dp.ps")

  # Run loop
  results <- list()

  for (window in initial_window:windows) {
    if (verbose) {
      message("Running window ", window, "/", windows, "...")
    }
    # Start and end of window
    start_nt <- 1 + stepsize * (window - 1)
    end_nt <- start_nt + winsize - 1

    # Allow for increased samples
    if (increased_sample) {
      start_nt <- ifelse(start_nt < 1, 1, start_nt)
      end_nt <- ifelse(end_nt > num_nt, num_nt, end_nt)
    }

    # Sequence contained in window
    win_seq <- seq[start_nt:end_nt] %>%
      .[!is.na(.)] %>%
      stringr::str_c(collapse = "")

    # Write temporary FASTA
    writeLines(
      text = c(paste(">", temp_seq_name), win_seq),
      con = temp_fasta_path
    )

    # Run RNAfold
    rnafold_res <- system(
      paste0(
        "cd ", temp_dir, "; ",
        "RNAfold ",
        rnafold_params,
        " -i ", temp_fasta_path
      ),
      intern = TRUE
    )

    # Parse Dot Plot Postcript file
    dotplot_data <- read_ps_dotplot(temp_dp_path)

    if (!is.null(dotplot_data)) {
      dotplot_data <- dotplot_data %>%
        dplyr::filter(box == "ubox") %>%
        dplyr::mutate(
          pos_i = pos_i + start_nt - 1,
          pos_j = pos_j + start_nt - 1,
          window = window
        ) %>%
        dplyr::select(window, pos_i, pos_j, prob)
      # dplyr::filter(abs(pos_j - pos_i) <= span)
    } else {
      next()
    }

    # Store results
    results[[as.character(window)]] <- list(
      "seq" = win_seq,
      "rnafold_res" = rnafold_res[-c(1, 2)],
      "dotplot_data" = dotplot_data
    )
    attr(results[[as.character(window)]], "start_nt") <- start_nt
    attr(results[[as.character(window)]], "end_nt") <- end_nt

    # Stop loop when end of sequence is reached
    if (!increased_sample & end_nt > num_nt) {
      break()
    }
  }

  # Additional results attributes
  attr(results, "seq_name") <- seq_name
  attr(results, "seq") <- stringr::str_c(seq, collapse = "")
  attr(results, "seq_length") <- num_nt

  return(results)
}
