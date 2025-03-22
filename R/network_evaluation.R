#' @title Calculate Network Prediction Performance Metrics
#'
#' @description
#' Calculates comprehensive performance metrics for evaluating predicted network structures,
#' including classification performance, precision-recall metrics, and network topology metrics.
#'
#' @md
#' @param network_table A data frame of predicted network structure containing:
#' \itemize{
#' \item *`regulator`* - Source nodes of the network edges
#' \item *`target`* - Target nodes of the network edges
#' \item *`weight`* - Edge weights representing prediction confidence
#' }
#' @param ground_truth A data frame of ground truth network with the same format as *`network_table`*.
#' @param metric_type The type of metric to return, default is *`all`*.
#' This can take any of the following choices:
#' \itemize{
#' \item *`all`* - Returns all available metrics with *Performance Metrics* plot
#' \item *`auc`* - Returns both AUROC and AUPRC with their plots
#' \item *`auroc`* - Area Under ROC Curve with plot
#' \item *`auprc`* - Area Under Precision-Recall Curve with plot
#' \item *`precision`* - Proportion of correct predictions among positive predictions
#' \item *`recall`* - Proportion of actual positives correctly identified
#' \item *`f1`* - Harmonic mean of precision and recall
#' \item *`accuracy`* - Overall classification accuracy
#' \item *`si`* - Set Intersection, counting correctly predicted edges
#' \item *`ji`* - Jaccard Index, measuring overlap between predicted and true networks
#' }
#' @param return_plot Logical value, default is *`FALSE`*, whether to generate visualization plots
#' @param line_color Color for plot lines, default is *`#1563cc`*
#' @param line_width Width for plot lines, default is *`1`*
#'
#' @return A list containing:
#' \itemize{
#' \item *`metrics`* - A data frame with requested metrics
#' \item *`plot`* - A plot object if return_plot = TRUE (optional)
#' }
#'
#' @export
#'
#' @examples
#' data("example_matrix")
#' data("example_ground_truth")
#' network_table <- inferCSN(example_matrix)
#' calculate_metrics(
#'   network_table,
#'   example_ground_truth,
#'   return_plot = TRUE
#' )
#' calculate_metrics(
#'   network_table,
#'   example_ground_truth,
#'   metric_type = "auroc"
#' )
calculate_metrics <- function(
    network_table,
    ground_truth,
    metric_type = c(
      "all",
      "auc",
      "auroc",
      "auprc",
      "precision",
      "recall",
      "f1",
      "accuracy",
      "si",
      "ji"
    ),
    return_plot = FALSE,
    line_color = "#1563cc",
    line_width = 1) {
  metric_type <- match.arg(metric_type)

  result <- switch(
    EXPR = metric_type,
    "all" = {
      auc_result <- calculate_auc(
        network_table,
        ground_truth,
        return_plot = FALSE
      )
      precision_result <- calculate_precision(
        network_table,
        ground_truth
      )
      recall_result <- calculate_recall(
        network_table,
        ground_truth
      )
      f1_result <- calculate_f1(
        network_table,
        ground_truth
      )
      accuracy_result <- calculate_accuracy(
        network_table,
        ground_truth
      )
      ji_result <- calculate_ji(
        network_table,
        ground_truth
      )
      si_result <- calculate_si(
        network_table,
        ground_truth
      )

      metrics_df <- purrr::map_dfr(
        list(
          auc_result,
          precision_result,
          recall_result,
          f1_result,
          accuracy_result,
          ji_result,
          si_result
        ), function(x) x$metrics
      )

      if (return_plot) {
        plot_metrics <- metrics_df[metrics_df$Metric != "SI", ]
        metrics_plot <- plot_all_metrics(plot_metrics, line_color)
        list(
          metrics = metrics_df,
          plot = metrics_plot
        )
      } else {
        list(
          metrics = metrics_df
        )
      }
    },
    "auc" = calculate_auc(
      network_table,
      ground_truth,
      return_plot,
      line_color,
      line_width
    ),
    "auroc" = calculate_auroc(
      network_table,
      ground_truth,
      return_plot,
      line_color,
      line_width
    ),
    "auprc" = calculate_auprc(
      network_table,
      ground_truth,
      return_plot,
      line_color,
      line_width
    ),
    "precision" = calculate_precision(
      network_table,
      ground_truth
    ),
    "recall" = calculate_recall(
      network_table,
      ground_truth
    ),
    "f1" = calculate_f1(
      network_table,
      ground_truth
    ),
    "accuracy" = calculate_accuracy(
      network_table,
      ground_truth
    ),
    "si" = calculate_si(
      network_table,
      ground_truth
    ),
    "ji" = calculate_ji(
      network_table,
      ground_truth
    )
  )

  return(result)
}

#' @title Prepare Binary Predictions
#'
#' @description
#' Prepares binary predictions and edge IDs from network tables
#'
#' @param network_table A data frame of predicted network structure
#' @param ground_truth A data frame of ground truth network
#'
#' @return A list containing binary predictions and edge information
#' @keywords internal
prepare_binary_predictions <- function(
    network_table,
    ground_truth) {
  gold <- prepare_calculate_metrics(
    network_table,
    ground_truth
  )

  results <- pROC::roc(
    gold$label ~ gold$weight,
    direction = "<",
    levels = c(0, 1)
  )
  select_value <- results$sensitivities + results$specificities - 1
  cut_value <- results$thresholds[select_value == max(select_value)][1]

  predictor_binary <- rep(0, length(results$predictor))
  predictor_binary[results$predictor >= cut_value] <- 1
  predictor_binary <- as.factor(predictor_binary) |>
    factor(levels = c(0, 1))
  true_label <- as.factor(gold$label) |>
    factor(levels = c(0, 1))

  pred_edges <- network_table[predictor_binary == 1, c("regulator", "target")]
  true_edges <- ground_truth[, c("regulator", "target")]
  pred_edge_ids <- paste(pred_edges$regulator, pred_edges$target, sep = "-")
  true_edge_ids <- paste(true_edges$regulator, true_edges$target, sep = "-")

  list(
    gold = gold,
    predictor_binary = predictor_binary,
    true_label = true_label,
    pred_edge_ids = pred_edge_ids,
    true_edge_ids = true_edge_ids
  )
}

.plot_auroc <- function(
    auc_curves,
    aurco,
    line_color,
    line_width) {
  subset(fortify(auc_curves), curvetype == "ROC") |>
    ggplot(aes(x = x, y = y)) +
    geom_line(color = line_color, linewidth = line_width) +
    geom_abline(
      slope = 1,
      intercept = 0,
      color = line_color,
      linetype = "dotted",
      linewidth = line_width
    ) +
    labs(
      title = paste("AUROC:", aurco),
      x = "False positive rate",
      y = "True positive rate"
    ) +
    theme_bw() +
    coord_fixed()
}

.plot_auprc <- function(
    auc_curves,
    auprc,
    line_color,
    line_width) {
  subset(fortify(auc_curves), curvetype == "PRC") |>
    ggplot(aes(x = x, y = y)) +
    geom_line(color = line_color, linewidth = line_width) +
    labs(
      title = paste("AUPRC:", auprc),
      x = "Recall",
      y = "Precision"
    ) +
    theme_bw() +
    coord_fixed()
}

#' @title Calculate AUROC Metric
#'
#' @description
#' Calculates AUROC metric with optional visualization
#'
#' @param network_table A data frame of predicted network structure
#' @param ground_truth A data frame of ground truth network
#' @param return_plot Logical value indicating whether to generate plot
#' @param line_color Color for plot lines
#' @param line_width Width for plot lines
#'
#' @return A list containing metric and optional plot
#' @export
#'
#' @examples
#' data("example_matrix")
#' data("example_ground_truth")
#' network_table <- inferCSN(example_matrix)
#' calculate_auroc(
#'   network_table,
#'   example_ground_truth,
#'   return_plot = TRUE
#' )
calculate_auroc <- function(
    network_table,
    ground_truth,
    return_plot = FALSE,
    line_color = "#1563cc",
    line_width = 1) {
  pred_data <- prepare_binary_predictions(
    network_table,
    ground_truth
  )
  gold <- pred_data$gold

  auc_curves <- precrec::evalmod(
    scores = gold$weight,
    labels = gold$label
  )
  auc <- attr(auc_curves, "auc")
  aurco <- as.numeric(sprintf("%0.3f", auc$aucs[1]))

  result <- data.frame(
    Metric = "AUROC",
    Value = aurco
  )

  if (return_plot) {
    auroc_plot <- .plot_auroc(
      auc_curves,
      aurco,
      line_color,
      line_width
    )

    return(
      list(
        metrics = result,
        plot = auroc_plot
      )
    )
  }

  return(
    list(
      metrics = result
    )
  )
}

#' @title Calculate AUPRC Metric
#'
#' @description
#' Calculates AUPRC metric with optional visualization
#'
#' @param network_table A data frame of predicted network structure
#' @param ground_truth A data frame of ground truth network
#' @param return_plot Logical value indicating whether to generate plot
#' @param line_color Color for plot lines
#' @param line_width Width for plot lines
#'
#' @return A list containing metric and optional plot
#' @export
#'
#' @examples
#' data("example_matrix")
#' data("example_ground_truth")
#' network_table <- inferCSN(example_matrix)
#' calculate_auprc(
#'   network_table,
#'   example_ground_truth,
#'   return_plot = TRUE
#' )
calculate_auprc <- function(
    network_table,
    ground_truth,
    return_plot = FALSE,
    line_color = "#1563cc",
    line_width = 1) {
  pred_data <- prepare_binary_predictions(
    network_table,
    ground_truth
  )
  gold <- pred_data$gold

  auc_curves <- precrec::evalmod(
    scores = gold$weight,
    labels = gold$label
  )
  auc <- attr(auc_curves, "auc")
  auprc <- as.numeric(sprintf("%0.3f", auc$aucs[2]))

  result <- data.frame(
    Metric = "AUPRC",
    Value = auprc
  )

  if (return_plot) {
    auprc_plot <- .plot_auprc(
      auc_curves,
      auprc,
      line_color,
      line_width
    )

    return(
      list(
        metrics = result,
        plot = auprc_plot
      )
    )
  }

  return(
    list(
      metrics = result
    )
  )
}

#' @title Calculate AUC Metrics
#'
#' @description
#' Calculates AUROC and AUPRC metrics with optional visualization
#'
#' @param network_table A data frame of predicted network structure
#' @param ground_truth A data frame of ground truth network
#' @param return_plot Logical value indicating whether to generate plots
#' @param line_color Color for plot lines
#' @param line_width Width for plot lines
#' @param tag_levels Tag levels for plot annotations
#'
#' @return A list containing metrics and optional plots
#' @export
#'
#' @examples
#' data("example_matrix")
#' data("example_ground_truth")
#' network_table <- inferCSN(example_matrix)
#' calculate_auc(
#'   network_table,
#'   example_ground_truth,
#'   return_plot = TRUE
#' )
calculate_auc <- function(
    network_table,
    ground_truth,
    return_plot = FALSE,
    line_color = "#1563cc",
    line_width = 1,
    tag_levels = "A") {
  auroc_result <- calculate_auroc(
    network_table,
    ground_truth,
    return_plot,
    line_color,
    line_width
  )
  auprc_result <- calculate_auprc(
    network_table,
    ground_truth,
    return_plot,
    line_color,
    line_width
  )

  result <- purrr::map_dfr(
    list(auroc_result, auprc_result),
    ~ .x$metrics
  )

  if (return_plot) {
    final_plot <- patchwork::wrap_plots(
      list(
        auroc_result$plot,
        auprc_result$plot
      ),
      ncol = 2
    ) + patchwork::plot_annotation(
      tag_levels = tag_levels
    )

    return(
      list(
        metrics = result,
        plot = final_plot
      )
    )
  }

  return(list(metrics = result))
}

#' @title Calculate Precision Metric
#'
#' @description
#' Calculates precision metric
#'
#' @param network_table A data frame of predicted network structure
#' @param ground_truth A data frame of ground truth network
#'
#' @return A list containing the metric
#' @export
#'
#' @examples
#' data("example_matrix")
#' data("example_ground_truth")
#' network_table <- inferCSN(example_matrix)
#' calculate_precision(
#'   network_table,
#'   example_ground_truth
#' )
calculate_precision <- function(
    network_table,
    ground_truth) {
  pred_data <- prepare_binary_predictions(
    network_table,
    ground_truth
  )

  pre <- as.vector(
    table(
      pred_data$predictor_binary,
      pred_data$true_label,
      dnn = c("Predicted", "Actual")
    )
  )

  TP <- pre[4]
  FP <- pre[2]

  precision <- as.numeric(sprintf("%0.3f", TP / (TP + FP)))

  list(
    metrics = data.frame(
      Metric = "Precision",
      Value = precision
    )
  )
}

#' @title Calculate Recall Metric
#'
#' @description
#' Calculates recall metric
#'
#' @param network_table A data frame of predicted network structure
#' @param ground_truth A data frame of ground truth network
#'
#' @return A list containing the metric
#' @export
#'
#' @examples
#' data("example_matrix")
#' data("example_ground_truth")
#' network_table <- inferCSN(example_matrix)
#' calculate_recall(
#'   network_table,
#'   example_ground_truth
#' )
calculate_recall <- function(
    network_table,
    ground_truth) {
  pred_data <- prepare_binary_predictions(
    network_table,
    ground_truth
  )

  pre <- as.vector(
    table(
      pred_data$predictor_binary,
      pred_data$true_label,
      dnn = c("Predicted", "Actual")
    )
  )

  TP <- pre[4]
  FN <- pre[3]

  recall <- as.numeric(sprintf("%0.3f", TP / (TP + FN)))

  list(
    metrics = data.frame(
      Metric = "Recall",
      Value = recall
    )
  )
}

#' @title Calculate F1 Score
#'
#' @description
#' Calculates F1 score
#'
#' @param network_table A data frame of predicted network structure
#' @param ground_truth A data frame of ground truth network
#'
#' @return A list containing the metric
#' @export
#'
#' @examples
#' data("example_matrix")
#' data("example_ground_truth")
#' network_table <- inferCSN(example_matrix)
#' calculate_f1(
#'   network_table,
#'   example_ground_truth
#' )
calculate_f1 <- function(
    network_table,
    ground_truth) {
  pred_data <- prepare_binary_predictions(
    network_table,
    ground_truth
  )

  pre <- as.vector(
    table(
      pred_data$predictor_binary,
      pred_data$true_label,
      dnn = c("Predicted", "Actual")
    )
  )

  TP <- pre[4]
  FP <- pre[2]
  FN <- pre[3]

  precision <- TP / (TP + FP)
  recall <- TP / (TP + FN)
  f1_score <- as.numeric(
    sprintf(
      "%0.3f",
      2 * (precision * recall) / (precision + recall)
    )
  )

  list(
    metrics = data.frame(
      Metric = "F1",
      Value = f1_score
    )
  )
}

#' @title Calculate Accuracy
#'
#' @description
#' Calculates accuracy metric
#'
#' @param network_table A data frame of predicted network structure
#' @param ground_truth A data frame of ground truth network
#'
#' @return A list containing the metric
#' @export
calculate_accuracy <- function(
    network_table,
    ground_truth) {
  pred_data <- prepare_binary_predictions(
    network_table,
    ground_truth
  )

  pre <- as.vector(
    table(
      pred_data$predictor_binary,
      pred_data$true_label,
      dnn = c("Predicted", "Actual")
    )
  )

  acc <- as.numeric(
    sprintf(
      "%0.3f",
      (pre[1] + pre[4]) / sum(pre)
    )
  )

  list(
    metrics = data.frame(
      Metric = "ACC",
      Value = acc
    )
  )
}

#' @title Calculate Set Intersection
#'
#' @description
#' Calculates Set Intersection (SI) metric
#'
#' @param network_table A data frame of predicted network structure
#' @param ground_truth A data frame of ground truth network
#'
#' @return A list containing the metric
#' @export
#'
#' @examples
#' data("example_matrix")
#' data("example_ground_truth")
#' network_table <- inferCSN(example_matrix)
#' calculate_si(
#'   network_table,
#'   example_ground_truth
#' )
calculate_si <- function(
    network_table,
    ground_truth) {
  pred_data <- prepare_binary_predictions(
    network_table,
    ground_truth
  )

  if (length(pred_data$pred_edge_ids) == 0 || length(pred_data$true_edge_ids) == 0) {
    si_score <- 0
  } else {
    si_score <- as.numeric(
      sprintf("%0.3f", length(intersect(pred_data$pred_edge_ids, pred_data$true_edge_ids)))
    )
  }

  list(
    metrics = data.frame(
      Metric = "SI",
      Value = si_score
    )
  )
}

#' @title Calculate Jaccard Index
#'
#' @description
#' Calculates Jaccard Index (JI) metric
#'
#' @param network_table A data frame of predicted network structure
#' @param ground_truth A data frame of ground truth network
#'
#' @return A list containing the metric
#' @export
#'
#' @examples
#' data("example_matrix")
#' data("example_ground_truth")
#' network_table <- inferCSN(example_matrix)
#' calculate_ji(
#'   network_table,
#'   example_ground_truth
#' )
calculate_ji <- function(
    network_table,
    ground_truth) {
  pred_data <- prepare_binary_predictions(
    network_table,
    ground_truth
  )

  if (length(pred_data$pred_edge_ids) == 0 || length(pred_data$true_edge_ids) == 0) {
    ji_score <- 0
  } else {
    intersection_size <- length(
      intersect(
        pred_data$pred_edge_ids,
        pred_data$true_edge_ids
      )
    )
    union_size <- length(
      union(
        pred_data$pred_edge_ids,
        pred_data$true_edge_ids
      )
    )
    ji_score <- as.numeric(
      sprintf(
        "%0.3f",
        intersection_size / union_size
      )
    )
  }

  list(
    metrics = data.frame(
      Metric = "JI",
      Value = ji_score
    )
  )
}

#' @title Create Performance Metrics Plot
#'
#' @description
#' Creates a bar plot of all performance metrics
#'
#' @param metrics_df Data frame containing metrics and their values
#' @param line_color Color for bars
#'
#' @return A ggplot object
#' @keywords internal
plot_all_metrics <- function(metrics_df, line_color = "#1563cc") {
  ggplot(
    metrics_df,
    aes(x = stats::reorder(Metric, -Value), y = Value)
  ) +
    geom_bar(stat = "identity", fill = line_color, width = 0.7) +
    geom_text(aes(label = Value), vjust = -0.5) +
    theme_bw() +
    labs(
      title = "Performance Metrics",
      x = "",
      y = "Score"
    ) +
    theme(
      axis.text.x = element_text(angle = 30, hjust = 1)
    ) +
    scale_y_continuous(
      expand = expansion(mult = c(0, 0.1))
    )
}

#' @title Network Edge Comparison Visualization
#'
#' @description
#' Generates visualizations comparing edges of two networks.
#'
#' @param network_table A data frame of predicted network structure.
#' @param ground_truth A data frame of ground truth network.
#' @param color_pattern A list of colors for different categories, with default values:
#' \itemize{
#'   \item *`predicted`* - Color for predicted edges (*`gray`*)
#'   \item *`ground_truth`* - Color for ground truth edges (*`#bb141a`*)
#'   \item *`overlap`* - Color for overlapping edges (*`#1966ad`*)
#'   \item *`total`* - Color for total counts (*`#6C757D`*)
#' }
#'
#' @return A patchwork plot object containing network edge comparison and distribution plots
#'
#' @export
#'
#' @examples
#' data("example_matrix")
#' data("example_ground_truth")
#' network_table <- inferCSN(example_matrix)
#' plot_edges_comparison(
#'   network_table,
#'   example_ground_truth
#' )
plot_edges_comparison <- function(
    network_table,
    ground_truth,
    color_pattern = list(
      predicted = "gray",
      ground_truth = "#bb141a",
      overlap = "#1966ad",
      total = "#6C757D"
    )) {
  all_genes <- sort(
    unique(
      c(
        network_table$regulator, network_table$target,
        ground_truth$regulator, ground_truth$target
      )
    )
  )

  pred_edges <- network_table[, c("regulator", "target")]
  true_edges <- ground_truth[, c("regulator", "target")]

  pred_edge_ids <- paste(
    pred_edges$regulator,
    pred_edges$target,
    sep = "-"
  )
  true_edge_ids <- paste(
    true_edges$regulator,
    true_edges$target,
    sep = "-"
  )
  overlap_edges <- intersect(
    pred_edge_ids,
    true_edge_ids
  )

  categories <- c(
    "Predicted Only" = "Predicted",
    "Ground Truth Only" = "Ground Truth",
    "Overlapping" = "Overlap"
  )

  edge_plot_data <- rbind(
    data.frame(
      regulator = pred_edges$regulator,
      target = pred_edges$target,
      type = "Predicted Only",
      stringsAsFactors = FALSE
    ),
    data.frame(
      regulator = true_edges$regulator,
      target = true_edges$target,
      type = "Ground Truth Only",
      stringsAsFactors = FALSE
    )
  )

  if (length(overlap_edges) > 0) {
    overlap_data <- data.frame(
      regulator = sub("-.*", "", overlap_edges),
      target = sub(".*-", "", overlap_edges),
      type = "Overlapping",
      stringsAsFactors = FALSE
    )
    edge_plot_data <- rbind(edge_plot_data, overlap_data)
  }

  edge_plot_data$type <- factor(
    edge_plot_data$type,
    levels = c(
      "Predicted Only",
      "Ground Truth Only",
      "Overlapping"
    )
  )

  edge_plot_data$regulator <- factor(
    edge_plot_data$regulator,
    levels = all_genes
  )
  edge_plot_data$target <- factor(
    edge_plot_data$target,
    levels = all_genes
  )

  network_plot <- ggplot() +
    geom_tile(
      data = edge_plot_data,
      aes(y = regulator, x = target, fill = type),
      color = "white",
      width = 0.9,
      height = 0.9
    ) +
    scale_fill_manual(
      values = c(
        "Predicted Only" = color_pattern$predicted,
        "Ground Truth Only" = color_pattern$ground_truth,
        "Overlapping" = color_pattern$overlap
      ),
      labels = categories
    ) +
    theme_bw() +
    labs(
      title = "Network Edge Comparison",
      y = "Regulator",
      x = "Target"
    ) +
    theme(
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        vjust = 1
      ),
      legend.position = "bottom",
      legend.title = element_blank(),
      panel.grid = element_line(color = "gray95"),
      plot.margin = margin(b = 20)
    ) +
    coord_fixed()

  intersection_size <- length(overlap_edges)
  predicted_only <- length(setdiff(pred_edge_ids, true_edge_ids))
  ground_truth_only <- length(setdiff(true_edge_ids, pred_edge_ids))
  total_predicted <- length(pred_edge_ids)
  total_ground_truth <- length(true_edge_ids)
  ordered_categories <- c(
    "Total Predicted",
    "Total Ground Truth",
    "Predicted Only",
    "Ground Truth Only",
    "Overlapping"
  )
  edge_details <- data.frame(
    Category = factor(
      ordered_categories,
      levels = ordered_categories
    ),
    count = c(
      total_predicted,
      total_ground_truth,
      predicted_only,
      ground_truth_only,
      intersection_size
    ),
    type = factor(
      c(
        "Total",
        "Total",
        "Predicted",
        "Ground Truth",
        "Overlap"
      ),
      levels = c(
        "Total",
        "Predicted",
        "Ground Truth",
        "Overlap"
      )
    )
  )

  edge_stats_plot <- ggplot(
    edge_details,
    aes(x = Category, y = count, fill = type)
  ) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_text(aes(label = count), vjust = -0.3) +
    scale_fill_manual(
      values = c(
        "Total" = color_pattern$total,
        "Predicted" = color_pattern$predicted,
        "Ground Truth" = color_pattern$ground_truth,
        "Overlap" = color_pattern$overlap
      )
    ) +
    theme_bw() +
    labs(
      title = "Edge Distribution",
      x = "",
      y = "Count"
    ) +
    theme(
      axis.text.x = element_text(angle = 30, hjust = 1),
      legend.position = "none",
      legend.title = element_blank()
    ) +
    scale_y_continuous(
      expand = expansion(mult = c(0, 0.1)),
      limits = function(x) c(0, max(x) * 1.05)
    )

  final_plot <- patchwork::wrap_plots(
    list(network_plot, edge_stats_plot),
    ncol = 2
  ) +
    patchwork::plot_annotation(tag_levels = "A")

  return(final_plot)
}
