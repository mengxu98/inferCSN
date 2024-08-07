% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/network_visualization.R
\name{plot_network_heatmap}
\alias{plot_network_heatmap}
\title{The heatmap of network}
\usage{
plot_network_heatmap(
  network_table,
  regulators = NULL,
  targets = NULL,
  switch_matrix = TRUE,
  show_names = FALSE,
  heatmap_size_lock = TRUE,
  heatmap_size = 5,
  heatmap_height = NULL,
  heatmap_width = NULL,
  heatmap_title = NULL,
  heatmap_color = c("#1966ad", "white", "#bb141a"),
  border_color = "gray",
  rect_color = NA,
  anno_width = 1,
  anno_height = 1,
  row_anno_type = NULL,
  column_anno_type = NULL,
  legend_name = "Weight",
  row_title = "Regulators"
)
}
\arguments{
\item{network_table}{The weight data table of network.}

\item{regulators}{Regulators list.}

\item{targets}{Targets list.}

\item{switch_matrix}{Logical value, default set to `TRUE`, whether to weight data table to matrix.}

\item{show_names}{Logical value, default set to `FALSE`, whether to show names of row and column.}

\item{heatmap_size_lock}{Lock the size of heatmap.}

\item{heatmap_size}{Default set to 5. The size of heatmap.}

\item{heatmap_height}{The height of heatmap.}

\item{heatmap_width}{The width of heatmap.}

\item{heatmap_title}{The title of heatmap.}

\item{heatmap_color}{Colors of heatmap.}

\item{border_color}{Default set to `gray`. Color of heatmap border.}

\item{rect_color}{Default set to `NA`. Color of heatmap rect.}

\item{anno_width}{Width of annotation.}

\item{anno_height}{Height of annotation.}

\item{row_anno_type}{Default set to `NULL`. c("boxplot", "barplot", "histogram", "density", "lines", "points", "horizon")}

\item{column_anno_type}{Default set to `NULL`. c("boxplot", "barplot", "histogram", "density", "lines", "points")}

\item{legend_name}{The name of legend.}

\item{row_title}{The title of row.}
}
\value{
Return a heatmap
}
\description{
The heatmap of network
}
\examples{
data("example_matrix")
data("example_ground_truth")
network_table <- inferCSN(example_matrix)

p1 <- plot_network_heatmap(
  example_ground_truth[, 1:3],
  heatmap_title = "Ground truth",
  legend_name = "Ground truth"
)
p2 <- plot_network_heatmap(
  network_table,
  heatmap_title = "inferCSN",
  legend_name = "inferCSN"
)
ComplexHeatmap::draw(p1 + p2)

p3 <- plot_network_heatmap(
  network_table,
  heatmap_title = "inferCSN",
  legend_name = "Weight1",
  heatmap_color = c("#20a485", "#410054", "#fee81f")
)
p4 <- plot_network_heatmap(
  network_table,
  heatmap_title = "inferCSN",
  legend_name = "Weight2",
  heatmap_color = c("#20a485", "white", "#fee81f")
)
ComplexHeatmap::draw(p3 + p4)

plot_network_heatmap(
  network_table,
  show_names = TRUE,
  rect_color = "gray90",
  row_anno_type = "density",
  column_anno_type = "barplot"
)

plot_network_heatmap(
  network_table,
  regulators = c("g1", "g2"),
  show_names = TRUE
)

plot_network_heatmap(
  network_table,
  targets = c("g1", "g2"),
  row_anno_type = "boxplot",
  column_anno_type = "histogram",
  show_names = TRUE
)

plot_network_heatmap(
  network_table,
  regulators = c("g1", "g3", "g5"),
  targets = c("g3", "g6", "g9"),
  show_names = TRUE
)
}
