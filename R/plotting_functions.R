#' @title Plot expression data
#'
#' @param data Input data.
#' @param smoothing_method Method for smoothing curve, `lm` or `loess`.
#' @param group_colors,title_color Group and title colors.
#' @param title,col_title,row_title,legend_title Plot titles.
#' @param legend_position The position of legend.
#' @param margins,marginal_type,margins_size Marginal plot controls.
#' @param compute_correlation,compute_correlation_method Correlation controls.
#' @param keep_aspect_ratio Whether to use a 1:1 aspect ratio.
#' @param facet Whether to facet by group.
#' @param se Whether to show smoothing uncertainty.
#' @param pointdensity Whether to show point density.
#' @return A ggplot object.
#' @export
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
  keep_aspect_ratio = TRUE,
  facet = FALSE,
  se = FALSE,
  pointdensity = TRUE
) {
  smoothing_method <- match.arg(
    smoothing_method,
    c("lm", "loess")
  )
  compute_correlation_method <- match.arg(
    compute_correlation_method,
    c("pearson", "spearman")
  )

  if (ncol(data) == 3) {
    colnames(data) <- c("x", "y", "cluster")
    p <- ggplot(
      data,
      aes(x = x, y = y, color = cluster)
    ) +
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
    p <- p +
      ggpubr::stat_cor(method = compute_correlation_method)
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
    marginal_type <- match.arg(
      marginal_type,
      c("density", "histogram", "boxplot", "violin", "densigram")
    )
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

#' @title Plot histogram
#'
#' @param data A numeric vector.
#' @param binwidth Width of the bins.
#' @param show_border Logical value, whether to show border of the bins.
#' @param border_color Color of the border.
#' @param alpha Alpha value of the bins.
#' @param theme Theme of the bins.
#' @param theme_begin Begin value of the theme.
#' @param theme_end End value of the theme.
#' @param theme_direction Direction of the theme.
#' @param legend_position The position of legend.
#'
#' @return A ggplot object
#' @export
#'
#' @examples
#' data(example_matrix)
#' network_table <- inferCSN(example_matrix)
#' plot_histogram(network_table[, 3])
plot_histogram <- function(
  data,
  binwidth = 0.01,
  show_border = FALSE,
  border_color = "black",
  alpha = 1,
  theme = "viridis",
  theme_begin = 0,
  theme_end = 0.5,
  theme_direction = -1,
  legend_position = "right"
) {
  data <- data.frame(weight = data)
  ggplot(data, aes(x = weight)) +
    geom_histogram(
      aes(fill = after_stat(count)),
      binwidth = binwidth,
      color = ifelse(show_border, border_color, NA),
      alpha = alpha
    ) +
    viridis::scale_fill_viridis(
      option = theme,
      begin = theme_begin,
      end = theme_end,
      direction = theme_direction
    ) +
    scale_x_continuous(name = "Weight") +
    scale_y_continuous(name = "Count") +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(
        color = "grey",
        size = 0.5
      ),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14, face = "bold"),
      legend.position = legend_position
    ) +
    theme_bw()
}

#' @title Plot an embedding
#'
#' @param matrix Input matrix.
#' @param labels Point labels.
#' @param method Embedding method.
#' @param colors Point colors.
#' @param point_size Point size.
#' @param seed Random seed.
#' @param cores Number of threads.
#' @return An embedding plot.
#' @export
plot_embedding <- function(
  matrix,
  labels = NULL,
  method = "pca",
  colors = RColorBrewer::brewer.pal(length(unique(labels)), "Set1"),
  seed = 1,
  point_size = 1,
  cores = 1
) {
  method <- match.arg(method, c("umap", "tsne", "pca"))

  set.seed(seed)
  result <- suppressMessages(
    switch(method,
      "umap" = {
        uwot::umap(
          matrix,
          n_components = 3,
          n_threads = cores,
          seed = seed
        )
      },
      "tsne" = {
        Rtsne::Rtsne(
          matrix,
          dims = 3,
          num_threads = cores
        )$Y
      },
      "pca" = {
        stats::prcomp(matrix, rank. = 3)$x
      }
    )
  ) |>
    as.data.frame()
  colnames(result) <- c("Dim1", "Dim2", "Dim3")

  if (!is.null(labels)) {
    result$label <- labels

    p <- plotly::plot_ly(
      data = result,
      x = ~Dim1,
      y = ~Dim2,
      z = ~Dim3,
      color = ~label,
      colors = colors,
      type = "scatter3d",
      mode = "markers",
      marker = list(size = point_size)
    )
  } else {
    p <- plotly::plot_ly(
      data = result,
      x = ~Dim1,
      y = ~Dim2,
      z = ~Dim3,
      type = "scatter3d",
      mode = "markers",
      marker = list(size = point_size)
    )
  }

  p <- plotly::layout(
    p,
    scene = list(
      xaxis = list(title = paste0(toupper(method), "_1")),
      yaxis = list(title = paste0(toupper(method), "_2")),
      zaxis = list(title = paste0(toupper(method), "_3"))
    ),
    title = paste0("3D ", toupper(method), " Visualization")
  )

  return(p)
}

#' @title Plot coefficients
#'
#' @md
#' @param data Input data.
#' @param style Plotting style: `"binary"`, `"gradient"`, or `"continuous"`.
#' @param positive_color Color for positive weights.
#' Default is `"#3d67a2"`.
#' @param negative_color Color for negative weights.
#' Default is `"#c82926"`.
#' @param neutral_color Color for weights near zero (used in `"continuous"` style).
#' Default is `"#cccccc"`.
#' @param bar_width Width of the bars.
#' Default is `0.7`.
#' @param text_size Size of the text for weight values.
#' Default is `3`.
#' @param show_values Whether to show weight values on bars.
#' Default is `TRUE`.
#'
#' @return A ggplot object
#' @export
#'
#' @examples
#' data(example_matrix)
#' network_table <- inferCSN(example_matrix, targets = "g1")
#' plot_coefficient(network_table)
#' plot_coefficient(network_table, style = "binary")
plot_coefficient <- function(
  data,
  style = "continuous",
  positive_color = "#3d67a2",
  negative_color = "#c82926",
  neutral_color = "#cccccc",
  bar_width = 0.7,
  text_size = 3,
  show_values = TRUE
) {
  p <- ggplot(
    data,
    aes(
      x = stats::reorder(regulator, weight),
      y = weight
    )
  ) +
    coord_flip() +
    labs(x = "Regulator", y = "Weight")

  if (style == "binary") {
    p <- p +
      geom_bar(stat = "identity", width = bar_width, aes(fill = weight > 0)) +
      scale_fill_manual(values = c(negative_color, positive_color))
  } else if (style == "continuous") {
    max_abs_weight <- max(abs(data$weight))
    p <- p +
      geom_bar(
        stat = "identity",
        width = bar_width,
        aes(fill = weight > 0, alpha = abs(weight) / max_abs_weight)
      ) +
      scale_fill_manual(values = c(negative_color, positive_color)) +
      scale_alpha_continuous(range = c(0.2, 1))
  } else {
    stop("Invalid style. Choose 'binary' or 'continuous'.")
  }

  p <- p +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "none",
      axis.text.y = element_text(size = 8),
      axis.text.x = element_text(size = 8),
      axis.title = element_text(size = 10, face = "bold")
    )

  if (show_values) {
    p <- p +
      geom_text(
        aes(
          label = sprintf("%.2f", weight),
          hjust = ifelse(weight > 0, -0.1, 1.1)
        ),
        size = text_size,
        color = "black"
      )
  }

  return(p)
}

#' @title Plot coefficients for multiple targets
#'
#' @param data Input data.
#' @param targets Targets to plot.
#' Default is `NULL`.
#' @param nrow Number of rows for the plot.
#' Default is `NULL`.
#' @param ... Other arguments passed to [plot_coefficient].
#'
#' @return A list of ggplot objects
#' @export
#'
#' @examples
#' data(example_matrix)
#' network_table <- inferCSN(
#'   example_matrix,
#'   targets = c("g1", "g2", "g3")
#' )
#' plot_coefficients(network_table, show_values = FALSE)
#' plot_coefficients(network_table, targets = "g1")
plot_coefficients <- function(
  data,
  targets = NULL,
  nrow = NULL,
  ...
) {
  if (is.null(targets)) {
    targets <- unique(data$target)
  }
  p_list <- lapply(targets, function(target) {
    plot_coefficient(data[data$target == target, ], ...)
  })
  if (!is.null(nrow)) {
    p <- patchwork::wrap_plots(p_list, nrow = nrow)
  } else {
    p <- patchwork::wrap_plots(p_list)
  }

  return(p)
}
