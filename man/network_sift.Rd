% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/network_sift.R
\name{network_sift}
\alias{network_sift}
\title{Sifting network}
\usage{
network_sift(
  network_table,
  matrix = NULL,
  meta_data = NULL,
  pseudotime_column = NULL,
  method = c("entropy", "max"),
  entropy_method = c("Shannon", "Renyi"),
  effective_entropy = FALSE,
  shuffles = 100,
  entropy_nboot = 300,
  lag_value = 1,
  entropy_p_value = 0.05,
  cores = 1,
  verbose = TRUE
)
}
\arguments{
\item{network_table}{The weight data table of network.}

\item{matrix}{The expression matrix.}

\item{meta_data}{The meta data for cells or samples.}

\item{pseudotime_column}{The column of pseudotime.}

\item{method}{The method used for filter edges.
Could be choose \code{entropy} or \code{max}.}

\item{entropy_method}{If setting \code{method} to \code{entropy},
could be choose \code{Shannon} or \code{Renyi} to compute entropy.}

\item{effective_entropy}{Default is \code{FALSE}.
Logical value, using effective entropy to filter weights or not.}

\item{shuffles}{Default is \code{100}.
The number of shuffles used to calculate the effective transfer entropy.}

\item{entropy_nboot}{Default is \code{300}.
The number of bootstrap replications for each direction of the estimated transfer entropy.}

\item{lag_value}{Default is \code{1}.
Markov order of gene expression values,
i.e. the number of lagged values affecting the current value of gene expression values.}

\item{entropy_p_value}{P value used to filter edges by entropy.}

\item{cores}{The number of cores to use for parallelization with \code{\link[foreach]{foreach}}, default is \emph{\code{1}}.}

\item{verbose}{Logical value, default is \emph{\code{TRUE}}, whether to print progress messages.}
}
\value{
Sifted network table
}
\description{
Sifting network
}
\examples{
\dontrun{
data("example_matrix")
data("example_meta_data")
data("example_ground_truth")

network_table <- inferCSN(example_matrix)
network_table_sifted <- network_sift(network_table)
network_table_sifted_entropy <- network_sift(
  network_table,
  matrix = example_matrix,
  meta_data = example_meta_data,
  pseudotime_column = "pseudotime",
  lag_value = 2,
  shuffles = 0,
  entropy_nboot = 0
)

plot_network_heatmap(
  example_ground_truth[, 1:3],
  heatmap_title = "Ground truth",
  show_names = TRUE,
  rect_color = "gray70"
)
plot_network_heatmap(
  network_table,
  heatmap_title = "Raw",
  show_names = TRUE,
  rect_color = "gray70"
)
plot_network_heatmap(
  network_table_sifted,
  heatmap_title = "Filtered",
  show_names = TRUE,
  rect_color = "gray70"
)
plot_network_heatmap(
  network_table_sifted_entropy,
  heatmap_title = "Filtered by entropy",
  show_names = TRUE,
  rect_color = "gray70"
)

calculate_metrics(
  network_table,
  example_ground_truth,
  return_plot = TRUE
)
calculate_metrics(
  network_table_sifted,
  example_ground_truth,
  return_plot = TRUE
)
calculate_metrics(
  network_table_sifted_entropy,
  example_ground_truth,
  return_plot = TRUE
)
}
}
