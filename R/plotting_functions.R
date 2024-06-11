#' @title plot_scatter
#'
#' @param data Input data
#' @param smoothing_method Method for smoothing curve, "lm" or "loess".
#' @param group_colors Colors for different groups.
#' @param title_color Color for the title.
#' @param title Main title for the plot.
#' @param col_title Title for the x-axis.
#' @param row_title Title for the y-axis.
#' @param legend_title Title for the legend.
#' @param legend_position The position of legend.
#' @param margins The position of marginal figure ("both", "x", "y").
#' @param marginal_type The type of marginal figure ("density", "histogram", "boxplot", "violin", "densigram").
#' @param margins_size The size of marginal figure, note the bigger size the smaller figure.
#' @param compute_correlation Whether to compute and print correlation on the figure.
#' @param compute_correlation_method Method to compute correlation ("pearson" or "spearman").
#' @param keep_aspect_ratio Logical value, whether to set aspect ratio to 1:1.
#' @param facet Faceting variable. If setting TRUE, all settings about margins will be inalidation.
#' @param se Display confidence interval around smooth.
#' @param pointdensity Plot point density when only provide 1 cluster.
#'
#' @return ggplot object
#' @export
#' @examples
#' data("example_matrix")
#' test_data <- data.frame(
#'   example_matrix[1:200, c(1, 7)],
#'   c = c(
#'     rep("c1", 40),
#'     rep("c2", 40),
#'     rep("c3", 40),
#'     rep("c4", 40),
#'     rep("c5", 40)
#'   )
#' )
#'
#' p1 <- plot_scatter(
#'   test_data,
#'   keep_aspect_ratio = TRUE
#' )
#' p2 <- plot_scatter(
#'   test_data,
#'   marginal_type = "boxplot",
#'   keep_aspect_ratio = TRUE
#' )
#' p1 + p2
#'
#' p3 <- plot_scatter(
#'   test_data,
#'   facet = TRUE,
#'   keep_aspect_ratio = TRUE
#' )
#' p3
#'
#' p4 <- plot_scatter(
#'   test_data[, 1:2],
#'   marginal_type = "histogram",
#'   keep_aspect_ratio = TRUE
#' )
#' p4
plot_scatter <- function(
    data,
    smoothing_method = "lm",
    group_colors = RColorBrewer::brewer.pal(9, "Set1"),
    title_color = "black",
    title = NULL,
    col_title = NULL,
    row_title = NULL,
    legend_title = NULL,
    legend_position = "bottom",
    margins = "both",
    marginal_type = NULL,
    margins_size = 10,
    compute_correlation = TRUE,
    compute_correlation_method = "pearson",
    keep_aspect_ratio = FALSE,
    facet = FALSE,
    se = FALSE,
    pointdensity = TRUE) {
  smoothing_method <- match.arg(smoothing_method, c("lm", "loess"))
  compute_correlation_method <- match.arg(
    compute_correlation_method, c("pearson", "spearman")
  )

  if (ncol(data) == 3) {
    colnames(data) <- c("x", "y", "cluster")
    p <- ggplot(data, aes(x = x, y = y, color = cluster)) +
      scale_color_manual(values = group_colors) +
      geom_point() +
      geom_smooth(
        method = smoothing_method,
        formula = "y ~ x",
        se = se
      )
    marginal_group_colour <- TRUE
    marginal_group_fill <- TRUE
  } else if (ncol(data) == 2) {
    colnames(data) <- c("x", "y")
    p <- ggplot(data, aes(x = x, y = y)) +
      geom_point(color = "#006699") +
      geom_smooth(
        method = smoothing_method,
        color = "#006699",
        formula = "y ~ x",
        se = se
      )
    if (pointdensity) {
      p <- p +
        ggpointdensity::geom_pointdensity() +
        viridis::scale_color_viridis()
    }
    marginal_group_colour <- FALSE
    marginal_group_fill <- FALSE
  } else {
    stop("Please provide a data frame with 2 or 3 columns.")
  }

  p <- p + theme_bw()

  if (compute_correlation) {
    p <- p + ggpubr::stat_cor(method = compute_correlation_method)
  }

  p <- p +
    labs(
      title = title,
      x = col_title,
      y = row_title,
      fill = legend_title
    ) +
    theme(
      plot.title = element_text(color = title_color),
      legend.position = legend_position
    )

  if (keep_aspect_ratio) {
    p <- p + coord_fixed(ratio = 1)
  }

  if (facet) {
    p <- p + facet_wrap(. ~ cluster)

    return(p)
  }

  if (!is.null(marginal_type)) {
    marginal_type <- match.arg(marginal_type, c("density", "histogram", "boxplot", "violin", "densigram"))
    margins <- match.arg(margins, c("both", "x", "y"))

    p <- suppressMessages(
      ggExtra::ggMarginal(
        p,
        margins = margins,
        type = marginal_type,
        groupColour = marginal_group_colour,
        groupFill = marginal_group_fill,
        xparams = list(colour = "#006699"),
        yparams = list(colour = "#006699"),
        size = margins_size
      )
    )
  }

  return(p)
}
