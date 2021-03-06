% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/base_pairs.R
\name{count_window_appearances}
\alias{count_window_appearances}
\title{Count Window Appearances}
\usage{
count_window_appearances(pos_i, pos_j, windows_list)
}
\arguments{
\item{pos_i}{Position of nucleotide \code{i}.}

\item{pos_j}{Position of nucleotide \code{j}.}

\item{windows_list}{List of \code{RNAfold} results for each sliding window. Result of \code{fold()} function.}
}
\value{
Integer of window appearance count.
}
\description{
Count on how many sliding windows the pair \code{(i, j)} was slid over.
}
