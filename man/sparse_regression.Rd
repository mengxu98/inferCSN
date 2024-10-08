% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sparse_regression.R
\name{sparse_regression}
\alias{sparse_regression}
\title{Sparse regression model}
\usage{
sparse_regression(
  x,
  y,
  cross_validation = FALSE,
  seed = 1,
  penalty = "L0",
  algorithm = "CD",
  regulators_num = ncol(x),
  n_folds = 10,
  subsampling = 1,
  r_threshold = 0,
  computation_method = "cor",
  verbose = TRUE,
  ...
)
}
\arguments{
\item{x}{The matrix of regulators.}

\item{y}{The vector of target.}

\item{cross_validation}{Logical value, default is \emph{\code{FALSE}}, whether to use cross-validation.}

\item{seed}{The random seed for cross-validation, default is \emph{\code{1}}.}

\item{penalty}{The type of regularization, default is \emph{\code{L0}}.
This can take either one of the following choices: \emph{\code{L0}}, \emph{\code{L0L1}}, and \emph{\code{L0L2}}.
For high-dimensional and sparse data, \emph{\code{L0L2}} is more effective.}

\item{algorithm}{The type of algorithm used to minimize the objective function, default is \emph{\code{CD}}.
Currently \emph{\code{CD}} and \emph{\code{CDPSI}} are supported.
The \emph{\code{CDPSI}} algorithm may yield better results, but it also increases running time.}

\item{regulators_num}{The number of non-zore coefficients, this value will affect the final performance.
The maximum support size at which to terminate the regularization path.
Recommend setting this to a small fraction of min(n,p) (e.g. 0.05 * min(n,p)) as L0 regularization typically selects a small portion of non-zeros.}

\item{n_folds}{The number of folds for cross-validation, default is \emph{\code{10}}.}

\item{subsampling}{The percent of all samples used for \code{\link{sparse_regression}}, default is \emph{\code{1}}.}

\item{r_threshold}{Threshold of \eqn{R^2} or correlation coefficient, default is \emph{\code{0}}.}

\item{computation_method}{The method used to compute `r``.}

\item{verbose}{Logical value, default is \emph{\code{TRUE}}, whether to print progress messages.}

\item{...}{Parameters for other methods.}
}
\value{
Coefficients
}
\description{
Sparse regression model
}
\examples{
data("example_matrix")
sparse_regression(
  example_matrix[, -1],
  example_matrix[, 1]
)
}
