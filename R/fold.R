#' Fold with Sliding Windows
#'
#' Perform folding inference using sliding windows and \code{RNAfold} from the
#' ViennaRNA Package 2.0 by Lorenz et al. <doi:10.1186/1748-7188-6-26>.
#'
#' @param filename The name of the file which the sequences in FASTA format are to be read from.
#' @param winsize Window size of the sliding windows.
#' @param stepsize Step size of the sliding windows.
#' @param same_num_samples If TRUE, fold will perform additional foldings at the beginning and end of the sequence. This allows the beginning and end of the sequence to have the same amount of samples/windows go through it. This is specially important for calculating Shannon's entropy.
#' @param rnafold_params Parameters used by \code{RNAfold}. The default is \code{'-p'}.
#' @param verbose If TRUE, fold will print information of performance.
#' @param rnafold_dir The path to the directory where RNAfold is. If NULL, it assumes that the dir is in the path
#'
#' @return List of \code{RNAfold} results for each sliding window.
#' @export
fold <- function(filename, winsize, stepsize, same_num_samples = TRUE, rnafold_params = "-p", verbose = FALSE, rnafold_dir = NULL) {
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
  if (same_num_samples) {
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
    if (same_num_samples) {
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
    if(!is.null(rnafold_dir)) { rnafold_cmd <- paste(rnafold_dir, "/RNAfold ", sep = "") }
    else { rnafold_cmd <- "RNAfold " }

    rnafold_res <- system(
      paste0(
        "cd ", temp_dir, "; ",
        rnafold_cmd,
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
    if (!same_num_samples & end_nt > num_nt) {
      break()
    }
  }

  # Additional results attributes
  attr(results, "seq_name") <- seq_name
  attr(results, "seq") <- stringr::str_c(seq, collapse = "")
  attr(results, "seq_length") <- num_nt

  return(results)
}
