% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/network_evaluation.R
\name{prepare_binary_predictions}
\alias{prepare_binary_predictions}
\title{Prepare Binary Predictions}
\usage{
prepare_binary_predictions(network_table, ground_truth)
}
\arguments{
\item{network_table}{A data frame of predicted network structure}

\item{ground_truth}{A data frame of ground truth network}
}
\value{
A list containing binary predictions and edge information
}
\description{
Prepares binary predictions and edge IDs from network tables
}
\keyword{internal}
