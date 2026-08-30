#' Plot network heatmaps
#'
#' @inheritParams network_format
#' @param switch_matrix Convert edge tables to matrices.
#' @param show_names,heatmap_size_lock Logical display controls.
#' @param show_names_position Name placement for multiple heatmaps.
#' @param heatmap_size,heatmap_height,heatmap_width Heatmap dimensions.
#' @param heatmap_title Heatmap title.
#' @param ncol,nrow Layout dimensions for multiple heatmaps.
#' @param performance_metrics Metrics appended to titles.
#' @param performance_ground_truth Ground-truth network for metrics.
#' @param truth_cell_border,truth_cell_border_legend Ground-truth border controls.
#' @param truth_match_border_color,truth_mismatch_border_color Border colors.
#' @param truth_cell_border_lwd Border width.
#' @param heatmap_palcolor,heatmap_palette Heatmap colors.
#' @param border_color,rect_color Border and cell colors.
#' @param row_anno,column_anno Enable annotations.
#' @param anno_width,anno_height Annotation dimensions.
#' @param row_anno_palette,row_anno_palcolor Row annotation colors.
#' @param col_anno_palette,col_anno_palcolor Column annotation colors.
#' @param row_anno_type,column_anno_type Annotation plot types.
#' @param legend_name,row_title Legend and row titles.
#' @return A heatmap object.
#' @export
plot_network_heatmap <- function(
  network_table,
  regulators = NULL,
  targets = NULL,
  switch_matrix = TRUE,
  show_names = FALSE,
  show_names_position = c("outer", "all"),
  rect_color = NA,
  heatmap_size_lock = TRUE,
  heatmap_size = 5,
  heatmap_height = NULL,
  heatmap_width = NULL,
  heatmap_title = NULL,
  ncol = NULL,
  nrow = NULL,
  performance_metrics = "none",
  performance_ground_truth = NULL,
  truth_cell_border = TRUE,
  truth_cell_border_legend = TRUE,
  truth_match_border_color = "#2E7D32",
  truth_mismatch_border_color = "#F2C94C",
  truth_cell_border_lwd = 1.2,
  border_color = "gray",
  heatmap_palette = "viridis",
  heatmap_palcolor = NULL,
  row_anno_palette = "Set1",
  row_anno_palcolor = NULL,
  col_anno_palette = "Set2",
  col_anno_palcolor = NULL,
  row_anno = FALSE,
  column_anno = FALSE,
  anno_width = 1,
  anno_height = 1,
  row_anno_type = c(
    "boxplot",
    "barplot",
    "histogram",
    "density",
    "lines",
    "points",
    "horizon"
  ),
  column_anno_type = c(
    "boxplot",
    "barplot",
    "histogram",
    "density",
    "lines",
    "points"
  ),
  legend_name = NULL,
  row_title = "Regulators"
) {
  show_names_position <- match.arg(show_names_position)

  if (
    is.list(network_table) &&
      !is.data.frame(network_table) &&
      !is.matrix(network_table)
  ) {
    network_tables <- network_table
    network_names <- names(network_tables)
    if (is.null(network_names) || any(!nzchar(network_names))) {
      network_names <- paste0("Network ", seq_along(network_tables))
    }
    ground_truth_index <- network_heatmap_ground_truth_index(network_names)

    performance_metrics <- network_heatmap_performance_metrics(
      performance_metrics
    )
    use_first_as_ground_truth <- length(performance_metrics) &&
      is.null(performance_ground_truth)
    if (use_first_as_ground_truth) {
      performance_ground_truth <- network_tables[[1]]
    }
    performance_title_labels <- network_heatmap_performance_title_labels(
      network_tables,
      ground_truth = performance_ground_truth,
      metrics = performance_metrics,
      skip = if (use_first_as_ground_truth) 1 else integer(0)
    )

    if (is.null(regulators)) {
      regulators <- valid_network_heatmap_names(unlist(
        lapply(network_tables, function(network_table) {
          rownames(network_to_heatmap_matrix(
            network_table,
            switch_matrix = switch_matrix
          ))
        }),
        use.names = FALSE
      ))
    }
    if (is.null(targets)) {
      targets <- valid_network_heatmap_names(unlist(
        lapply(network_tables, function(network_table) {
          colnames(network_to_heatmap_matrix(
            network_table,
            switch_matrix = switch_matrix
          ))
        }),
        use.names = FALSE
      ))
    }

    weight_matrices <- lapply(
      network_tables,
      network_to_aligned_heatmap_matrix,
      regulators = regulators,
      targets = targets,
      switch_matrix = switch_matrix
    )
    names(weight_matrices) <- network_names
    truth_border_matrix <- NULL
    truth_border_summaries <- NULL
    if (isTRUE(truth_cell_border)) {
      truth_table <- performance_ground_truth
      if (is.null(truth_table) && !is.na(ground_truth_index)) {
        truth_table <- network_tables[[ground_truth_index]]
      }
      if (!is.null(truth_table)) {
        truth_border_matrix <- network_to_aligned_heatmap_matrix(
          truth_table,
          regulators = regulators,
          targets = targets,
          switch_matrix = switch_matrix
        )
        truth_border_summaries <- lapply(seq_along(weight_matrices), function(i) {
          network_heatmap_truth_cell_summary(
            weight_matrices[[i]],
            truth_matrix = truth_border_matrix,
            skip = !is.na(ground_truth_index) && i == ground_truth_index
          )
        })
      }
    }
    name_visibility <- network_heatmap_name_visibility(
      n = length(weight_matrices),
      ncol = ncol,
      nrow = nrow,
      position = show_names_position
    )

    legend_names <- vapply(
      seq_along(network_tables),
      function(i) {
        if (is.null(legend_name)) {
          network_value_label(
            network_tables[[i]],
            switch_matrix = switch_matrix
          )
        } else {
          select_plot_label(legend_name, i, NULL) %||%
            network_value_label(
              network_tables[[i]],
              switch_matrix = switch_matrix
            )
        }
      },
      character(1)
    )
    heatmap_names <- make.unique(legend_names, sep = "_")

    heatmaps <- lapply(seq_along(weight_matrices), function(i) {
      color_function <- network_heatmap_col_fun(
        as.vector(weight_matrices[[i]]),
        heatmap_palcolor = heatmap_palcolor,
        heatmap_palette = heatmap_palette
      )
      build_network_heatmap(
        weight_matrix = weight_matrices[[i]],
        color_function = color_function,
        cell_border_matrix = network_heatmap_truth_cell_borders(
          weight_matrices[[i]],
          truth_matrix = truth_border_matrix,
          skip = !is.na(ground_truth_index) && i == ground_truth_index,
          match_color = truth_match_border_color,
          mismatch_color = truth_mismatch_border_color
        ),
        cell_border_lwd = truth_cell_border_lwd,
        show_names = show_names,
        show_row_names = show_names && name_visibility$row[[i]],
        show_column_names = show_names && name_visibility$column[[i]],
        heatmap_size_lock = heatmap_size_lock,
        heatmap_size = heatmap_size,
        heatmap_height = heatmap_height,
        heatmap_width = heatmap_width,
        heatmap_title = append_network_heatmap_performance_label(
          select_plot_label(heatmap_title, i, network_names[[i]]),
          performance_title_labels[[i]]
        ),
        border_color = border_color,
        rect_color = rect_color,
        anno_width = anno_width,
        anno_height = anno_height,
        row_anno_palette = row_anno_palette,
        row_anno_palcolor = row_anno_palcolor,
        col_anno_palette = col_anno_palette,
        col_anno_palcolor = col_anno_palcolor,
        row_anno = row_anno,
        column_anno = column_anno,
        row_anno_type = row_anno_type,
        column_anno_type = column_anno_type,
        heatmap_name = heatmap_names[[i]],
        legend_name = legend_names[[i]],
        row_title = row_title
      )
    })

    extra_legends <- lapply(seq_along(heatmaps), function(i) {
      show_truth_legend <- isTRUE(truth_cell_border_legend)
      legend_summary <- if (is.null(truth_border_summaries)) {
        NULL
      } else {
        truth_border_summaries[[i]]
      }
      invisible_placeholder <- FALSE
      if (
        !is.na(ground_truth_index) &&
          i == ground_truth_index &&
          !is.null(truth_border_summaries)
      ) {
        non_truth <- setdiff(seq_along(truth_border_summaries), ground_truth_index)
        non_empty <- non_truth[vapply(
          truth_border_summaries[non_truth],
          function(x) sum(x) > 0,
          logical(1)
        )]
        if (length(non_empty)) {
          legend_summary <- truth_border_summaries[[non_empty[[1]]]]
          invisible_placeholder <- TRUE
          show_truth_legend <- TRUE
        }
      }
      truth_border_legend <- network_heatmap_truth_cell_legend(
        legend_summary,
        show = show_truth_legend,
        match_color = truth_match_border_color,
        mismatch_color = truth_mismatch_border_color,
        visible = !invisible_placeholder
      )
      if (is.null(truth_border_legend)) {
        NULL
      } else {
        list(truth_border_legend)
      }
    })

    p <- arrange_network_heatmaps(
      heatmaps,
      ncol = ncol,
      nrow = nrow,
      extra_legends = extra_legends
    )
    return(p)
  }

  performance_metrics <- network_heatmap_performance_metrics(
    performance_metrics
  )
  performance_title_label <- ""
  if (length(performance_metrics)) {
    if (is.null(performance_ground_truth)) {
      stop(
        "performance_ground_truth must be provided when performance_metrics ",
        "is not \"none\" for a single network_table.",
        call. = FALSE
      )
    }
    performance_title_label <- network_heatmap_performance_label(
      network_table,
      ground_truth = performance_ground_truth,
      metrics = performance_metrics
    )
  }

  weight_matrix <- network_to_aligned_heatmap_matrix(
    network_table,
    regulators = regulators,
    targets = targets,
    switch_matrix = switch_matrix
  )

  color_function <- network_heatmap_col_fun(
    as.vector(weight_matrix),
    heatmap_palcolor = heatmap_palcolor,
    heatmap_palette = heatmap_palette
  )
  legend_name <- legend_name %||%
    network_value_label(network_table, switch_matrix = switch_matrix)
  single_truth_matrix <- if (!is.null(performance_ground_truth)) {
    network_to_aligned_heatmap_matrix(
      performance_ground_truth,
      regulators = rownames(weight_matrix),
      targets = colnames(weight_matrix),
      switch_matrix = switch_matrix
    )
  } else {
    NULL
  }
  single_truth_summary <- list(network_heatmap_truth_cell_summary(
    weight_matrix,
    truth_matrix = single_truth_matrix,
    skip = !isTRUE(truth_cell_border)
  ))

  p <- build_network_heatmap(
    weight_matrix = weight_matrix,
    color_function = color_function,
    cell_border_matrix = network_heatmap_truth_cell_borders(
      weight_matrix,
      truth_matrix = single_truth_matrix,
      skip = !isTRUE(truth_cell_border),
      match_color = truth_match_border_color,
      mismatch_color = truth_mismatch_border_color
    ),
    cell_border_lwd = truth_cell_border_lwd,
    show_names = show_names,
    show_row_names = show_names,
    show_column_names = show_names,
    heatmap_size_lock = heatmap_size_lock,
    heatmap_size = heatmap_size,
    heatmap_height = heatmap_height,
    heatmap_width = heatmap_width,
    heatmap_title = append_network_heatmap_performance_label(
      heatmap_title,
      performance_title_label
    ),
    border_color = border_color,
    rect_color = rect_color,
    anno_width = anno_width,
    anno_height = anno_height,
    row_anno_palette = row_anno_palette,
    row_anno_palcolor = row_anno_palcolor,
    col_anno_palette = col_anno_palette,
    col_anno_palcolor = col_anno_palcolor,
    row_anno = row_anno,
    column_anno = column_anno,
    row_anno_type = row_anno_type,
    column_anno_type = column_anno_type,
    heatmap_name = legend_name,
    legend_name = legend_name,
    row_title = row_title
  )

  truth_legend <- network_heatmap_truth_cell_legend(
    single_truth_summary,
    show = isTRUE(truth_cell_border_legend),
    match_color = truth_match_border_color,
    mismatch_color = truth_mismatch_border_color
  )
  if (!is.null(truth_legend)) {
    p <- arrange_network_heatmaps(
      list(p),
      ncol = 1,
      extra_legends = list(list(truth_legend))
    )
  }

  return(p)
}

network_to_heatmap_matrix <- function(
  network_table,
  regulators = NULL,
  targets = NULL,
  switch_matrix = TRUE
) {
  if (switch_matrix) {
    network_table <- as.data.frame(network_table, stringsAsFactors = FALSE)
    if (ncol(network_table) < 3) {
      stop("network_table must contain at least three columns.")
    }
    network_table <- network_table[, 1:3, drop = FALSE]
    colnames(network_table) <- c("row", "col", "value")
    valid_edges <- valid_network_heatmap_name(network_table$row) &
      valid_network_heatmap_name(network_table$col)
    network_table <- network_table[valid_edges, , drop = FALSE]
    if (!nrow(network_table)) {
      stop("network_table contains no valid regulator-target pairs.")
    }
    if (!is.numeric(network_table$value)) {
      value <- as.character(network_table$value)
      mapped <- ifelse(
        value %in% c("+", "Activation", "activation", "activating", "positive"),
        1,
        ifelse(
          value %in%
            c("-", "Repression", "repression", "repressing", "negative"),
          -1,
          suppressWarnings(as.numeric(value))
        )
      )
      mapped[is.na(mapped)] <- 0
      network_table$value <- mapped
    }
    weight_matrix <- thisutils::table_to_matrix(
      network_table,
      regulators,
      targets
    )
  } else {
    weight_matrix <- as.matrix(network_table)
    if (!is.null(rownames(weight_matrix))) {
      valid_rows <- valid_network_heatmap_name(rownames(weight_matrix))
      weight_matrix <- weight_matrix[valid_rows, , drop = FALSE]
    }
    if (!is.null(colnames(weight_matrix))) {
      valid_cols <- valid_network_heatmap_name(colnames(weight_matrix))
      weight_matrix <- weight_matrix[, valid_cols, drop = FALSE]
    }
  }
  weight_matrix <- filter_sort_matrix(
    weight_matrix,
    regulators = regulators,
    targets = targets
  )
  storage.mode(weight_matrix) <- "numeric"
  weight_matrix
}

network_to_aligned_heatmap_matrix <- function(
  network_table,
  regulators = NULL,
  targets = NULL,
  switch_matrix = TRUE
) {
  weight_matrix <- network_to_heatmap_matrix(
    network_table,
    regulators = NULL,
    targets = NULL,
    switch_matrix = switch_matrix
  )
  if (is.null(regulators)) {
    regulators <- rownames(weight_matrix)
  }
  if (is.null(targets)) {
    targets <- colnames(weight_matrix)
  }
  regulators <- valid_network_heatmap_names(regulators)
  targets <- valid_network_heatmap_names(targets)
  if (!length(regulators) || !length(targets)) {
    stop("network_table contains no valid regulator-target pairs.")
  }
  aligned_matrix <- matrix(
    0,
    nrow = length(regulators),
    ncol = length(targets),
    dimnames = list(regulators, targets)
  )
  row_index <- intersect(regulators, rownames(weight_matrix))
  col_index <- intersect(targets, colnames(weight_matrix))
  aligned_matrix[row_index, col_index] <- weight_matrix[row_index, col_index]
  aligned_matrix
}

valid_network_heatmap_name <- function(x) {
  !is.na(x) & nzchar(as.character(x))
}

valid_network_heatmap_names <- function(x) {
  x <- as.character(x)
  x <- x[valid_network_heatmap_name(x)]
  x[!duplicated(x)]
}

network_heatmap_performance_metrics <- function(metrics) {
  if (is.null(metrics) || !length(metrics)) {
    return(character(0))
  }
  metrics <- unique(tolower(as.character(metrics)))
  metrics <- metrics[nzchar(metrics)]
  if (!length(metrics) || identical(metrics, "none")) {
    return(character(0))
  }
  metrics <- setdiff(metrics, "none")
  allowed <- c(
    "all",
    "all_no_epr",
    "auc",
    "auroc",
    "auprc",
    "epr",
    "precision",
    "recall",
    "f1",
    "accuracy",
    "si",
    "ji"
  )
  invalid <- setdiff(metrics, allowed)
  if (length(invalid)) {
    stop(
      "performance_metrics contains unsupported metric(s): ",
      paste(invalid, collapse = ", "),
      call. = FALSE
    )
  }
  metrics
}

network_heatmap_performance_title_labels <- function(
  network_tables,
  ground_truth,
  metrics,
  skip = integer(0)
) {
  labels <- rep("", length(network_tables))
  if (!length(metrics)) {
    return(labels)
  }
  if (is.null(ground_truth)) {
    stop("performance_ground_truth must be provided.", call. = FALSE)
  }
  for (i in seq_along(network_tables)) {
    if (i %in% skip) {
      next
    }
    labels[[i]] <- network_heatmap_performance_label(
      network_tables[[i]],
      ground_truth = ground_truth,
      metrics = metrics
    )
  }
  labels
}

network_heatmap_performance_label <- function(
  network_table,
  ground_truth,
  metrics
) {
  metric_results <- lapply(metrics, function(metric) {
    calculate_metrics(
      network_table = network_table,
      ground_truth = ground_truth,
      metric_type = metric,
      return_plot = FALSE
    )$metrics
  })
  metric_results <- do.call(rbind, metric_results)
  metric_results <- metric_results[!duplicated(metric_results$Metric), ]
  metric_labels <- paste0(
    metric_results$Metric,
    ": ",
    format_network_heatmap_metric_value(metric_results$Value)
  )
  paste(metric_labels, collapse = ", ")
}

format_network_heatmap_metric_value <- function(value, digits = 3) {
  value <- as.numeric(value)
  ifelse(
    is.na(value),
    "NA",
    format(round(value, digits), trim = TRUE, scientific = FALSE)
  )
}

append_network_heatmap_performance_label <- function(title, label) {
  if (is.null(label) || !nzchar(label)) {
    return(title)
  }
  if (is.null(title) || is.na(title) || !nzchar(title)) {
    return(paste0("(", label, ")"))
  }
  paste0(title, " (", label, ")")
}

network_heatmap_legend_breaks <- function(weight_matrix) {
  values <- as.numeric(weight_matrix)
  values <- values[is.finite(values)]
  if (!length(values)) {
    return(0)
  }
  min_weight <- min(values)
  max_weight <- max(values)
  if (min_weight < 0 && max_weight > 0) {
    max_abs <- max(abs(c(min_weight, max_weight)))
    return(unique(c(-max_abs, 0, max_abs)))
  }
  if (identical(min_weight, max_weight)) {
    return(min_weight)
  }
  unique(c(min_weight, mean(c(min_weight, max_weight)), max_weight))
}

format_network_heatmap_legend_labels <- function(values, digits = 3) {
  values <- as.numeric(values)
  labels <- format(signif(values, digits), trim = TRUE, scientific = FALSE)
  labels[labels == "-0"] <- "0"
  labels
}

network_heatmap_col_fun <- function(
  values,
  heatmap_palcolor,
  heatmap_palette = NULL
) {
  values <- values[is.finite(values)]
  if (!length(values)) {
    values <- 0
  }
  min_weight <- min(values)
  max_weight <- max(values)
  if (identical(min_weight, max_weight)) {
    min_weight <- min_weight - 1
    max_weight <- max_weight + 1
  }
  pal_colors <- network_palette_colors(
    palcolor = heatmap_palcolor,
    palette = heatmap_palette,
    n = 101,
    type = "continuous"
  )
  if (min_weight >= 0) {
    color_function <- circlize::colorRamp2(
      c(min_weight, max_weight),
      pal_colors[c(ceiling(length(pal_colors) / 2), length(pal_colors))]
    )
  } else if (max_weight <= 0) {
    color_function <- circlize::colorRamp2(
      c(min_weight, max_weight),
      pal_colors[c(1, ceiling(length(pal_colors) / 2))]
    )
  } else {
    max_abs <- max(abs(c(min_weight, max_weight)))
    color_function <- circlize::colorRamp2(
      c(-max_abs, 0, max_abs),
      pal_colors[c(1, ceiling(length(pal_colors) / 2), length(pal_colors))]
    )
  }
  color_function
}

network_palette_colors <- function(
  palcolor = NULL,
  palette = NULL,
  n = 101,
  type = c("continuous", "discrete")
) {
  type <- match.arg(type)
  n <- max(1, as.integer(n))
  if (!is.null(palcolor)) {
    if (type == "continuous") {
      n <- max(3, n)
    }
    colors <- grDevices::colorRampPalette(palcolor)(n)
  } else if (is.null(palette)) {
    colors <- thisplot::palette_colors(
      x = seq_len(n),
      n = n,
      type = type
    )
  } else {
    colors <- thisplot::palette_colors(
      x = seq_len(n),
      n = n,
      palette = palette,
      type = type
    )
  }
  colors <- unname(colors)
  if (type == "continuous" && length(colors) < 3) {
    colors <- grDevices::colorRampPalette(colors)(3)
  }
  colors
}

select_plot_label <- function(label, index, default) {
  if (is.null(label)) {
    return(default)
  }
  if (length(label) >= index) {
    return(label[[index]])
  }
  label[[1]]
}

arrange_network_heatmaps <- function(
  heatmaps,
  ncol = NULL,
  nrow = NULL,
  extra_grob = NULL,
  extra_legends = NULL
) {
  if (
    is.null(ncol) &&
      is.null(nrow) &&
      is.null(extra_grob) &&
      is.null(extra_legends)
  ) {
    return(Reduce(`+`, heatmaps))
  }
  if (is.null(extra_legends)) {
    extra_legends <- rep(list(NULL), length(heatmaps))
  }

  layout <- if (is.null(ncol) && is.null(nrow)) {
    list(ncol = length(heatmaps), nrow = 1)
  } else {
    network_heatmap_grid_layout(
      n = length(heatmaps),
      ncol = ncol,
      nrow = nrow
    )
  }

  grobs <- lapply(seq_along(heatmaps), function(i) {
    grid::grid.grabExpr(
      ComplexHeatmap::draw(
        heatmaps[[i]],
        newpage = FALSE,
        heatmap_legend_list = extra_legends[[i]]
      )
    )
  })

  children <- grid::gList()
  for (i in seq_along(grobs)) {
    row_index <- ((i - 1) %/% layout$ncol) + 1
    col_index <- ((i - 1) %% layout$ncol) + 1
    children <- grid::gList(
      children,
      grid::editGrob(
        grobs[[i]],
        vp = grid::viewport(
          layout.pos.row = row_index,
          layout.pos.col = col_index
        )
      )
    )
  }
  if (!is.null(extra_grob)) {
    children <- grid::gList(
      children,
      grid::grobTree(
        extra_grob,
        vp = grid::viewport(
          layout.pos.row = seq_len(layout$nrow),
          layout.pos.col = layout$ncol + 1
        )
      )
    )
  }

  layout_widths <- if (is.null(extra_grob)) {
    grid::unit(rep(1, layout$ncol), "null")
  } else {
    grid::unit.c(
      grid::unit(rep(1, layout$ncol), "null"),
      grid::unit(3.4, "cm")
    )
  }
  p <- grid::gTree(
    children = children,
    vp = grid::viewport(
      layout = grid::grid.layout(
        nrow = layout$nrow,
        ncol = layout$ncol + as.integer(!is.null(extra_grob)),
        widths = layout_widths
      )
    ),
    name = "infercsn_heatmap_grid"
  )
  class(p) <- c("infercsn_heatmap_grid", class(p))
  p
}

network_heatmap_grid_layout <- function(n, ncol = NULL, nrow = NULL) {
  if (!is.null(ncol)) {
    ncol <- as.integer(ncol)
    if (length(ncol) != 1 || is.na(ncol) || ncol < 1) {
      stop("ncol must be a positive integer.", call. = FALSE)
    }
  }
  if (!is.null(nrow)) {
    nrow <- as.integer(nrow)
    if (length(nrow) != 1 || is.na(nrow) || nrow < 1) {
      stop("nrow must be a positive integer.", call. = FALSE)
    }
  }

  if (is.null(ncol)) {
    ncol <- ceiling(n / nrow)
  }
  if (is.null(nrow)) {
    nrow <- ceiling(n / ncol)
  }
  if (ncol * nrow < n) {
    stop("ncol * nrow must be at least the number of heatmaps.", call. = FALSE)
  }

  list(ncol = ncol, nrow = nrow)
}

network_heatmap_name_visibility <- function(
  n,
  ncol = NULL,
  nrow = NULL,
  position = c("outer", "all")
) {
  position <- match.arg(position)
  if (position == "all" || n <= 1) {
    return(list(
      row = rep(TRUE, n),
      column = rep(TRUE, n)
    ))
  }

  layout <- if (is.null(ncol) && is.null(nrow)) {
    list(ncol = n, nrow = 1)
  } else {
    network_heatmap_grid_layout(n = n, ncol = ncol, nrow = nrow)
  }
  row_index <- ((seq_len(n) - 1) %/% layout$ncol) + 1
  col_index <- ((seq_len(n) - 1) %% layout$ncol) + 1
  rightmost_col <- vapply(
    row_index,
    function(row) {
      max(col_index[row_index == row])
    },
    numeric(1)
  )
  bottom_row <- vapply(
    col_index,
    function(col) {
      max(row_index[col_index == col])
    },
    numeric(1)
  )

  list(
    row = col_index == rightmost_col,
    column = row_index == bottom_row
  )
}

network_heatmap_ground_truth_index <- function(network_names) {
  if (is.null(network_names)) {
    return(NA_integer_)
  }
  hit <- grep(
    "^(ground[ _-]*truth|truth)$",
    network_names,
    ignore.case = TRUE
  )
  if (length(hit)) {
    return(hit[[1]])
  }
  NA_integer_
}

network_heatmap_truth_cell_borders <- function(
  weight_matrix,
  truth_matrix = NULL,
  skip = FALSE,
  match_color = "#2E7D32",
  mismatch_color = "#F2C94C"
) {
  if (isTRUE(skip) || is.null(truth_matrix)) {
    return(NULL)
  }
  truth_matrix <- truth_matrix[
    rownames(weight_matrix),
    colnames(weight_matrix),
    drop = FALSE
  ]
  pred_sign <- sign(weight_matrix)
  truth_sign <- sign(truth_matrix)
  hit <- truth_sign != 0 & pred_sign != 0
  border_matrix <- matrix(
    NA_character_,
    nrow = nrow(weight_matrix),
    ncol = ncol(weight_matrix),
    dimnames = dimnames(weight_matrix)
  )
  border_matrix[hit & pred_sign == truth_sign] <- match_color
  border_matrix[hit & pred_sign != truth_sign] <- mismatch_color
  border_matrix
}

network_heatmap_truth_cell_summary <- function(
  weight_matrix,
  truth_matrix = NULL,
  skip = FALSE
) {
  empty <- c(
    truth_edges = 0,
    predicted_edges = 0,
    hits = 0,
    true_direction = 0,
    false_direction = 0
  )
  if (isTRUE(skip) || is.null(truth_matrix)) {
    return(empty)
  }
  truth_matrix <- truth_matrix[
    rownames(weight_matrix),
    colnames(weight_matrix),
    drop = FALSE
  ]
  pred_sign <- sign(weight_matrix)
  truth_sign <- sign(truth_matrix)
  hit <- truth_sign != 0 & pred_sign != 0
  c(
    truth_edges = sum(truth_sign != 0, na.rm = TRUE),
    predicted_edges = sum(pred_sign != 0, na.rm = TRUE),
    hits = sum(hit, na.rm = TRUE),
    true_direction = sum(hit & pred_sign == truth_sign, na.rm = TRUE),
    false_direction = sum(hit & pred_sign != truth_sign, na.rm = TRUE)
  )
}

network_heatmap_truth_cell_legend <- function(
  summaries,
  show = TRUE,
  match_color = "#2E7D32",
  mismatch_color = "#F2C94C",
  visible = TRUE
) {
  if (!isTRUE(show) || is.null(summaries) || !length(summaries)) {
    return(NULL)
  }
  if (is.numeric(summaries) && !is.null(names(summaries))) {
    summaries <- list(summaries)
  }
  summary <- Reduce(`+`, summaries)
  if (is.null(summary) || !sum(summary)) {
    return(NULL)
  }
  summary[["truth_edges"]] <- max(vapply(
    summaries,
    function(x) x[["truth_edges"]],
    numeric(1)
  ))

  ComplexHeatmap::Legend(
    title = if (isTRUE(visible)) "Truth border" else "Truth border",
    labels = c(
      paste0("Truth: ", summary[["truth_edges"]]),
      paste0("Pred: ", summary[["predicted_edges"]]),
      paste0("Hits: ", summary[["hits"]]),
      paste0("TRUE: ", summary[["true_direction"]]),
      paste0("FALSE: ", summary[["false_direction"]])
    ),
    type = "points",
    pch = if (isTRUE(visible)) c(15, 15, 15, 0, 0) else rep(15, 5),
    legend_gp = grid::gpar(
      col = if (isTRUE(visible)) {
        c("white", "white", "white", match_color, mismatch_color)
      } else {
        rep("#FFFFFF00", 5)
      },
      fill = if (isTRUE(visible)) {
        c("white", "white", "white", NA, NA)
      } else {
        rep("#FFFFFF00", 5)
      },
      lwd = 1.2
    ),
    title_gp = grid::gpar(
      fontface = "plain",
      fontsize = 8,
      col = if (isTRUE(visible)) "black" else "#FFFFFF00"
    ),
    labels_gp = grid::gpar(
      fontsize = 7,
      col = if (isTRUE(visible)) "black" else "#FFFFFF00"
    ),
    background = if (isTRUE(visible)) "white" else "#FFFFFF00",
    ncol = 1
  )
}

network_heatmap_cell_border_fun <- function(border_matrix = NULL, lwd = 1.2) {
  if (is.null(border_matrix)) {
    return(NULL)
  }
  force(border_matrix)
  force(lwd)
  function(j, i, x, y, width, height, fill) {
    border_color <- border_matrix[i, j]
    if (!is.na(border_color) && nzchar(border_color)) {
      grid::grid.rect(
        x = x,
        y = y,
        width = width,
        height = height,
        gp = grid::gpar(fill = NA, col = border_color, lwd = lwd)
      )
    }
  }
}

#' @export
print.infercsn_heatmap_grid <- function(x, ...) {
  grid::grid.newpage()
  grid::grid.draw(x)
  invisible(x)
}

register_infercsn_heatmap_grid_draw_method <- function() {
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
    return(invisible(FALSE))
  }

  methods::setOldClass(c("infercsn_heatmap_grid", "gTree", "grob", "gDesc"))
  methods::setMethod(
    ComplexHeatmap::draw,
    methods::signature(object = "infercsn_heatmap_grid"),
    function(object, newpage = TRUE, ...) {
      if (isTRUE(newpage)) {
        grid::grid.newpage()
      }
      grid::grid.draw(object)
      invisible(object)
    }
  )
  invisible(TRUE)
}

.onLoad <- function(libname, pkgname) {
  register_infercsn_heatmap_grid_draw_method()
}

network_value_label <- function(network_table, switch_matrix = TRUE) {
  if (!switch_matrix) {
    return("value")
  }

  network_table <- as.data.frame(network_table, stringsAsFactors = FALSE)
  value_name <- colnames(network_table)[[3]] %||% ""
  if (nzchar(value_name)) {
    return(value_name)
  }

  "value"
}

build_network_heatmap <- function(
  weight_matrix,
  color_function,
  cell_border_matrix = NULL,
  cell_border_lwd = 1.2,
  show_names = FALSE,
  show_row_names = show_names,
  show_column_names = show_names,
  heatmap_size_lock = TRUE,
  heatmap_size = 5,
  heatmap_height = NULL,
  heatmap_width = NULL,
  heatmap_title = NULL,
  border_color = "gray",
  rect_color = NA,
  anno_width = 1,
  anno_height = 1,
  row_anno_palette = NULL,
  row_anno_palcolor = NULL,
  col_anno_palette = NULL,
  col_anno_palcolor = NULL,
  row_anno = FALSE,
  column_anno = FALSE,
  row_anno_type = c(
    "boxplot",
    "barplot",
    "histogram",
    "density",
    "lines",
    "points",
    "horizon"
  ),
  column_anno_type = c(
    "boxplot",
    "barplot",
    "histogram",
    "density",
    "lines",
    "points"
  ),
  heatmap_name = legend_name,
  legend_name = NULL,
  row_title = "Regulators"
) {
  unique_regulators <- rownames(weight_matrix)
  unique_targets <- colnames(weight_matrix)

  if (show_names) {
    if (is.null(heatmap_height) || is.null(heatmap_width)) {
      heatmap_height <- length(unique_regulators) / 2
      heatmap_width <- length(unique_targets) / 2
    }
  } else {
    if (is.null(heatmap_height) || is.null(heatmap_width)) {
      heatmap_height <- heatmap_size *
        length(unique_regulators) /
        length(unique_targets)
      heatmap_width <- heatmap_size
    }
  }

  if (isTRUE(row_anno)) {
    row_anno_fill <- network_palette_colors(
      palcolor = row_anno_palcolor,
      palette = row_anno_palette,
      n = nrow(weight_matrix),
      type = "discrete"
    )
    col_anno_fill <- network_palette_colors(
      palcolor = col_anno_palcolor,
      palette = col_anno_palette,
      n = ncol(weight_matrix),
      type = "discrete"
    )
    row_anno_type <- match.arg(row_anno_type)
    row_anno <- switch(row_anno_type,
      "boxplot" = ComplexHeatmap::rowAnnotation(
        " " = ComplexHeatmap::anno_boxplot(
          weight_matrix,
          width = grid::unit(anno_width, "cm"),
          gp = grid::gpar(fill = row_anno_fill)
        )
      ),
      "barplot" = ComplexHeatmap::rowAnnotation(
        " " = ComplexHeatmap::anno_barplot(
          abs(weight_matrix),
          width = grid::unit(anno_width, "cm"),
          gp = grid::gpar(fill = col_anno_fill)
        )
      ),
      "histogram" = ComplexHeatmap::rowAnnotation(
        " " = ComplexHeatmap::anno_histogram(
          weight_matrix,
          width = grid::unit(anno_width, "cm"),
          gp = grid::gpar(fill = row_anno_fill)
        )
      ),
      "density" = ComplexHeatmap::rowAnnotation(
        " " = ComplexHeatmap::anno_density(
          weight_matrix,
          width = grid::unit(anno_width, "cm"),
          gp = grid::gpar(fill = row_anno_fill)
        )
      ),
      "lines" = ComplexHeatmap::rowAnnotation(
        " " = ComplexHeatmap::anno_lines(
          weight_matrix,
          width = grid::unit(anno_width, "cm"),
          gp = grid::gpar(fill = row_anno_fill)
        )
      ),
      "points" = ComplexHeatmap::rowAnnotation(
        " " = ComplexHeatmap::anno_points(
          weight_matrix,
          width = grid::unit(anno_width, "cm"),
          gp = grid::gpar(fill = row_anno_fill)
        )
      ),
      "horizon" = ComplexHeatmap::rowAnnotation(
        " " = ComplexHeatmap::anno_horizon(
          weight_matrix,
          width = grid::unit(anno_width, "cm"),
          gp = grid::gpar(fill = row_anno_fill)
        )
      )
    )
  } else {
    row_anno <- NULL
  }

  if (isTRUE(column_anno)) {
    row_anno_fill <- network_palette_colors(
      palcolor = row_anno_palcolor,
      palette = row_anno_palette,
      n = nrow(weight_matrix),
      type = "discrete"
    )
    col_anno_fill <- network_palette_colors(
      palcolor = col_anno_palcolor,
      palette = col_anno_palette,
      n = ncol(weight_matrix),
      type = "discrete"
    )
    column_anno_type <- match.arg(column_anno_type)
    column_anno <- switch(column_anno_type,
      "boxplot" = ComplexHeatmap::columnAnnotation(
        " " = ComplexHeatmap::anno_boxplot(
          weight_matrix,
          height = grid::unit(anno_height, "cm"),
          gp = grid::gpar(fill = col_anno_fill)
        )
      ),
      "barplot" = ComplexHeatmap::columnAnnotation(
        " " = ComplexHeatmap::anno_barplot(
          abs(weight_matrix),
          height = grid::unit(anno_height, "cm"),
          gp = grid::gpar(fill = row_anno_fill)
        )
      ),
      "histogram" = ComplexHeatmap::columnAnnotation(
        " " = ComplexHeatmap::anno_histogram(
          weight_matrix,
          height = grid::unit(anno_height, "cm"),
          gp = grid::gpar(fill = col_anno_fill)
        )
      ),
      "density" = ComplexHeatmap::columnAnnotation(
        " " = ComplexHeatmap::anno_density(
          weight_matrix,
          height = grid::unit(anno_height, "cm"),
          gp = grid::gpar(fill = col_anno_fill)
        )
      ),
      "lines" = ComplexHeatmap::columnAnnotation(
        " " = ComplexHeatmap::anno_lines(
          weight_matrix,
          height = grid::unit(anno_height, "cm"),
          gp = grid::gpar(fill = col_anno_fill)
        )
      ),
      "points" = ComplexHeatmap::columnAnnotation(
        " " = ComplexHeatmap::anno_points(
          weight_matrix,
          height = grid::unit(anno_height, "cm"),
          gp = grid::gpar(fill = col_anno_fill)
        )
      )
    )
  } else {
    column_anno <- NULL
  }

  if (heatmap_size_lock) {
    width <- grid::unit(heatmap_width, "cm")
    height <- grid::unit(heatmap_height, "cm")
  } else {
    width <- NULL
    height <- NULL
  }

  legend_breaks <- network_heatmap_legend_breaks(weight_matrix)
  cell_fun <- network_heatmap_cell_border_fun(
    cell_border_matrix,
    lwd = cell_border_lwd
  )
  p <- ComplexHeatmap::Heatmap(
    weight_matrix,
    name = heatmap_name,
    col = color_function,
    heatmap_legend_param = list(
      title = legend_name,
      title_gp = grid::gpar(fontface = "plain", fontsize = 7),
      labels_gp = grid::gpar(fontface = "plain", fontsize = 7),
      legend_height = grid::unit(17.8, "mm"),
      grid_width = grid::unit(2.7, "mm"),
      at = legend_breaks,
      labels = format_network_heatmap_legend_labels(legend_breaks),
      border = NA
    ),
    column_title = heatmap_title,
    row_title = row_title,
    column_title_gp = grid::gpar(fontsize = 11),
    row_title_gp = grid::gpar(fontsize = 11),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = show_row_names,
    show_column_names = show_column_names,
    row_names_gp = grid::gpar(fontface = "italic", fontsize = 10),
    column_names_gp = grid::gpar(fontface = "italic", fontsize = 10),
    row_names_max_width = grid::unit(1.25, "cm"),
    column_names_rot = 45,
    gap = grid::unit(1.8, "mm"),
    border = border_color,
    rect_gp = grid::gpar(col = rect_color),
    width = width,
    height = height,
    top_annotation = column_anno,
    left_annotation = row_anno,
    cell_fun = cell_fun
  )

  return(p)
}

#' @title Plot dynamic networks
#'
#' @inheritParams network_format
#' @param legend_position The position of legend.
#'
#' @return A ggplot2 object
#' @export
#'
#' @examples
#' data(example_matrix)
#' network_table <- inferCSN(example_matrix)
#' example_edge <- network_table[1, ]
#' plot_static_networks(
#'   network_table,
#'   regulators = example_edge$regulator
#' )
#' plot_static_networks(
#'   network_table,
#'   targets = example_edge$target
#' )
#' plot_static_networks(
#'   network_table,
#'   regulators = example_edge$regulator,
#'   targets = example_edge$target
#' )
plot_static_networks <- function(
  network_table,
  regulators = NULL,
  targets = NULL,
  legend_position = "right"
) {
  network_table <- network_format(
    network_table,
    regulators = regulators,
    targets = targets
  )

  net <- igraph::graph_from_data_frame(
    network_table[, c("regulator", "target", "weight", "Interaction")],
    directed = FALSE
  )

  layout <- igraph::layout_with_fr(net)
  rownames(layout) <- igraph::V(net)$name
  layout_ordered <- layout[igraph::V(net)$name, ]
  regulator_network <- ggnetwork(
    net,
    layout = layout_ordered,
    cell.jitter = 0
  )

  regulator_network$is_regulator <- as.character(
    regulator_network$name %in% regulators
  )
  cols <- c("Activation" = "#3366cc", "Repression" = "#ff0066")

  g <- ggplot() +
    geom_edges(
      data = regulator_network,
      aes(
        x = x,
        y = y,
        xend = xend,
        yend = yend,
        linewidth = weight,
        color = Interaction
      ),
      curvature = 0.1,
      alpha = .6
    ) +
    geom_nodes(
      data = regulator_network[regulator_network$is_regulator == "FALSE", ],
      aes(x = x, y = y),
      color = "darkgray",
      size = 3,
      alpha = .5
    ) +
    geom_nodes(
      data = regulator_network[regulator_network$is_regulator == "TRUE", ],
      aes(x = x, y = y),
      color = "#8C4985",
      size = 6,
      alpha = .8
    ) +
    scale_color_manual(values = cols) +
    geom_nodelabel_repel(
      data = regulator_network[regulator_network$is_regulator == "FALSE", ],
      aes(x = x, y = y, label = name),
      size = 2,
      color = "#5A8BAD"
    ) +
    geom_nodelabel_repel(
      data = regulator_network[regulator_network$is_regulator == "TRUE", ],
      aes(x = x, y = y, label = name),
      size = 3.5,
      color = "black"
    ) +
    theme_blank() +
    theme(legend.position = legend_position)

  return(g)
}

#' @title Plot contrast networks
#'
#' @md
#' @inheritParams plot_static_networks
#' @param degree_value Degree value to filter nodes.
#' Default is `0`.
#' @param weight_value Weight value to filter edges.
#' Default is `0`.
#' @param legend_position The position of legend.
#' Default is `"bottom"`.
#'
#' @return
#' A ggplot2 object.
#' @export
#'
#' @examples
#' data(example_matrix)
#' network_table <- inferCSN(example_matrix)
#' plot_contrast_networks(utils::head(network_table, 50))
plot_contrast_networks <- function(
  network_table,
  degree_value = 0,
  weight_value = 0,
  legend_position = "bottom"
) {
  graph <- network_format(network_table) |>
    tidygraph::as_tbl_graph() |>
    dplyr::mutate(
      degree = tidygraph::centrality_degree(mode = "out")
    ) |>
    dplyr::filter(degree > degree_value) |>
    tidygraph::activate(edges)

  g <- ggraph(graph, layout = "linear", circular = TRUE) +
    geom_edge_arc(
      aes(
        colour = Interaction,
        filter = weight > weight_value,
        edge_width = weight
      ),
      arrow = arrow(length = unit(3, "mm")),
      start_cap = square(3, "mm"),
      end_cap = circle(3, "mm")
    ) +
    scale_edge_width(range = c(0, 1)) +
    facet_edges(~Interaction) +
    geom_node_point(aes(size = degree), colour = "#A1B7CE") +
    geom_node_text(aes(label = name), repel = TRUE) +
    coord_fixed() +
    theme_graph(
      base_family = "serif",
      foreground = "steelblue",
      fg_text_colour = "white"
    ) +
    theme(legend.position = legend_position)

  return(g)
}

#' Plot dynamic networks
#'
#' @inheritParams network_format
#' @param celltypes_order Cell-type order.
#' @param ntop Number of top genes to plot.
#' @param title Figure title.
#' @param theme_type Plot theme.
#' @param plot_type Output type.
#' @param layout Network layout.
#' @param nrow Number of rows.
#' @param figure_save Whether to save the figure.
#' @param figure_name Output filename.
#' @param figure_width,figure_height Figure dimensions.
#' @param seed Random seed.
#' @return A dynamic network plot.
#' @export
plot_dynamic_networks <- function(
  network_table,
  celltypes_order,
  ntop = 10,
  title = NULL,
  theme_type = "theme_void",
  plot_type = "ggplot",
  layout = "fruchtermanreingold",
  nrow = 2,
  figure_save = FALSE,
  figure_name = NULL,
  figure_width = 6,
  figure_height = 6,
  seed = 1
) {
  names(network_table) <- c("regulator", "target", "weight", "celltype")
  network_table$regulator <- as.character(network_table$regulator)
  network_table$target <- as.character(network_table$target)
  network_table <- network_table[
    network_table$regulator != network_table$target,
  ]
  celltypes_list <- unique(intersect(celltypes_order, network_table$celltype))

  network_table <- purrr::map_dfr(
    celltypes_list,
    .f = function(x) {
      network_table[which(network_table$celltype == x), ]
    }
  )

  nodes <- unique(
    c(network_table$regulator, network_table$target)
  )
  dnodes <- data.frame(id = 1:length(nodes), label = nodes)
  edges <- dplyr::left_join(
    network_table,
    dnodes,
    by = c("regulator" = "label")
  ) |>
    dplyr::rename(from = id) |>
    dplyr::left_join(
      dnodes,
      by = c("target" = "label")
    ) |>
    dplyr::rename(to = id) |>
    dplyr::select(from, to, weight, celltype)
  edges$Interaction <- ifelse(
    edges$weight > 0,
    "Activation",
    "Repression"
  )
  edges$weight <- abs(edges$weight)

  dedges <- unique(edges)
  dnodes$label <- gsub("\\.", "-", dnodes$label)
  network_data <- network::network(
    dedges,
    vertex.attr = dnodes,
    matrix.type = "edgelist",
    ignore.eval = FALSE,
    directed = TRUE,
    multiple = TRUE
  )

  set.seed(seed)
  layout <- match.arg(
    layout,
    c("fruchtermanreingold", "kamadakawai")
  )
  ggnetwork_data <- ggnetwork(
    network_data,
    arrow.size = 0.1,
    arrow.gap = 0.015,
    by = "celltype",
    weights = "weight",
    layout = layout
  )

  nodes_data <- purrr::map_dfr(
    celltypes_list,
    .f = function(x) {
      nodes_data_celltype <- network_table[
        which(network_table$celltype == x),
      ] |>
        dplyr::group_by(
          regulator
        ) |>
        dplyr::summarise(
          targets_num = dplyr::n()
        ) |>
        dplyr::arrange(
          dplyr::desc(targets_num)
        ) |>
        as.data.frame()
      nodes_data_celltype$label_genes <- as.character(
        nodes_data_celltype$regulator
      )
      if (nrow(nodes_data_celltype) > ntop) {
        cf <- nodes_data_celltype$targets_num[ntop]
        nodes_data_celltype$label_genes[which(
          nodes_data_celltype$targets_num < cf
        )] <- ""
      } else if (nrow(nodes_data_celltype) == 0) {
        return()
      }
      nodes_data_celltype$celltype <- x

      return(nodes_data_celltype)
    }
  )

  names(nodes_data)[1] <- "label"
  ggnetwork_data <- merge(
    ggnetwork_data,
    nodes_data,
    by = c("label", "celltype"),
    all.x = T
  )
  ggnetwork_data$targets_num[which(is.na(ggnetwork_data$targets_num))] <- 0
  ggnetwork_data$label_genes[which(is.na(ggnetwork_data$label_genes))] <- ""
  ggnetwork_data$celltype <- factor(
    ggnetwork_data$celltype,
    levels = celltypes_order
  )
  cols <- c("Activation" = "#3366cc", "Repression" = "#ff0066")

  plot_type <- match.arg(
    plot_type,
    c("ggplot", "animate", "ggplotly")
  )
  p <- ggplot(ggnetwork_data, aes(x, y, xend = xend, yend = yend))
  if (plot_type == "ggplotly") {
    p <- p +
      geom_edges(
        aes(color = Interaction),
        linewidth = 0.7,
        arrow = arrow(length = unit(3, "pt"), type = "closed")
      )
  } else {
    p <- p +
      geom_edges(
        aes(color = Interaction, alpha = weight),
        linewidth = 0.7,
        arrow = arrow(length = unit(3, "pt"), type = "closed")
      )
  }
  p <- p +
    geom_nodes(
      aes(size = targets_num),
      color = "darkgray",
      alpha = 0.9
    ) +
    geom_nodetext(
      aes(label = label_genes),
      color = "black"
    ) +
    theme(aspect.ratio = 2, legend.position = "bottom") +
    scale_color_manual(values = cols)

  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }

  theme_type <- match.arg(
    theme_type,
    c("theme_void", "theme_facet", "theme_blank")
  )
  p <- switch(
    EXPR = theme_type,
    "theme_void" = p + theme_void(),
    "theme_blank" = p + theme_blank(),
    "theme_facet" = p + theme_facet()
  )
  if (plot_type == "ggplot") {
    p <- p +
      facet_wrap(~celltype, nrow = nrow)
    if (figure_save) {
      if (is.null(figure_name)) {
        figure_name <- "networks.pdf"
      }
      ggsave(
        figure_name,
        p,
        width = figure_width,
        height = figure_height
      )
    }
  }

  if (plot_type == "animate") {
    p <- p + gganimate::transition_states(states = celltype)
    p <- gganimate::animate(
      p,
      render = gganimate::gifski_renderer()
    )
    if (figure_save) {
      if (is.null(figure_name)) {
        figure_name <- "networks.gif"
      }
      gganimate::anim_save(figure_name, animation = p)
    }
  }

  if (plot_type == "ggplotly") {
    p <- p +
      facet_wrap(~celltype, nrow = nrow)
    p <- plotly::ggplotly(p)
  }

  return(p)
}
