% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/subsampling.R
\name{subsampling}
\alias{subsampling}
\title{Subsampling function}
\usage{
subsampling(
  matrix,
  subsampling_method = c("sample", "meta_cells", "pseudobulk"),
  subsampling_ratio = 1,
  seed = 1,
  verbose = TRUE,
  ...
)
}
\arguments{
\item{matrix}{The input matrix to be subsampled.}

\item{subsampling_method}{The method to use for subsampling. Options are "sample", "pseudobulk" or "meta_cells".}

\item{subsampling_ratio}{The percent of all samples used for \code{\link{fit_srm}}, default is \emph{\code{1}}.}

\item{seed}{The random seed for cross-validation, default is \emph{\code{1}}.}

\item{verbose}{Logical value, default is \emph{\code{TRUE}}, whether to print progress messages.}

\item{...}{Parameters for other methods.}
}
\value{
The subsampled matrix.
}
\description{
This function subsamples a matrix using either random sampling or meta cells method.
}
\examples{
data("example_matrix")
data("example_ground_truth")
subsample_matrix <- subsampling(
  example_matrix,
  subsampling_ratio = 0.5
)
subsample_matrix_2 <- subsampling(
  example_matrix,
  subsampling_method = "meta_cells",
  subsampling_ratio = 0.5,
  fast_pca = FALSE
)
subsample_matrix_3 <- subsampling(
  example_matrix,
  subsampling_method = "pseudobulk",
  subsampling_ratio = 0.5
)

calculate_metrics(
  inferCSN(example_matrix),
  example_ground_truth,
  return_plot = TRUE
)
calculate_metrics(
  inferCSN(subsample_matrix),
  example_ground_truth,
  return_plot = TRUE
)
calculate_metrics(
  inferCSN(subsample_matrix_2),
  example_ground_truth,
  return_plot = TRUE
)
calculate_metrics(
  inferCSN(subsample_matrix_3),
  example_ground_truth,
  return_plot = TRUE
)
}
