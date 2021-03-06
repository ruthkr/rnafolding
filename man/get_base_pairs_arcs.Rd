% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/base_pairs.R
\name{get_base_pairs_arcs}
\alias{get_base_pairs_arcs}
\title{Calculate Base-Pairs Arc Trajectories}
\usage{
get_base_pairs_arcs(base_pairs_data)
}
\arguments{
\item{base_pairs_data}{Data frame specifying folding probability of the base-pair \code{(pos_i, pos_j)}.}
}
\value{
Data frame arc trajectories for each base-pair.
}
\description{
Calculate base-pairs segments to build base-pairing probabilities plot.
}
