#' @title Calculate Network Prediction Performance Metrics
#'
#' @description
#' Calculates comprehensive performance metrics for evaluating predicted network structures,
#' including classification performance, precision-recall metrics, and network topology metrics.
#'
#' @md
#' @param network_table A data frame of predicted network structure containing:
#'
#' * `regulator` - Source nodes of the network edges
#'
#' * `target` - Target nodes of the network edges
#'
#' * `weight` - Edge weights representing prediction confidence
#'
#' @param ground_truth A data frame of ground truth network with the same format as *`network_table`*.
#'
#' @param type The type of metric to return, default is *`all`*.
#' This can take any of the following choices:
#'
#' * *`all`* - Returns all available metrics
#'
#' * *`auroc`* - Area Under ROC Curve, evaluating overall classification ability
#'
#' * *`auprc`* - Area Under Precision-Recall Curve, focusing on positive class prediction
#'
#' * *`precision`* - Proportion of correct predictions among positive predictions
#'
#' * *`recall`* - Proportion of actual positives correctly identified
#'
#' * *`f1`* - Harmonic mean of precision and recall
#'
#' * *`si`* - Set Intersection, counting correctly predicted edges
#'
#' * *`ji`* - Jaccard Index, measuring overlap between predicted and true networks
#'
#' @param plot Logical value, default is *`FALSE`*, whether to generate visualization plots.
#' If *`TRUE`*, generates:
#'
#' * ROC curve plot for AUROC evaluation
#'
#' * Precision-Recall curve plot for AUPRC evaluation
#'
#' * Confusion matrix heatmap for classification results
#'
#' * Network comparison plot showing edge overlap
#'
#' * Metrics summary bar plot
#'
#' @param line_color The color for plot lines, default is *`#1563cc`*.
#' @param line_width The width for plot lines, default is *`1`*.
#'
#' @return
#' * If type="all": A data frame with all metrics
#' * If specific type: A single numeric value with the requested metric
#' * If plot=TRUE: Displays relevant visualization plots
#'
#' @examples
#' data("example_matrix")
#' data("example_ground_truth")
#' inferCSN(example_matrix) |>
#'   calculate_metrics(example_ground_truth, plot = TRUE)
#'
#' @export
calculate_metrics <- function(
    network_table,
    ground_truth,
    type = "all",
    plot = FALSE,
    line_color = "#1563cc",
    line_width = 1) {
  gold <- prepare_calculate_metrics(
    network_table,
    ground_truth
  )

  results <- pROC::roc(
    gold$label ~ gold$weight,
    direction = "<",
    levels = c(0, 1)
  )

  auc_curves <- precrec::evalmod(
    scores = gold$weight,
    labels = gold$label
  )
  auc <- attr(auc_curves, "auc")
  aurco <- as.numeric(
    sprintf("%0.3f", auc$aucs[1])
  )
  auprc <- as.numeric(
    sprintf("%0.3f", auc$aucs[2])
  )

  select_value <- results$sensitivities + results$specificities - 1
  cut_value <- results$thresholds[select_value == max(select_value)][1]

  predictor_binary <- rep(0, length(results$predictor))
  predictor_binary[results$predictor >= cut_value] <- 1
  predictor_binary <- as.factor(predictor_binary) |>
    factor(levels = c(0, 1))
  reverse_label <- as.factor(gold$label) |>
    factor(levels = c(0, 1))

  pre <- as.vector(
    table(
      predictor_binary,
      reverse_label
    )
  )
  acc <- as.numeric(
    sprintf("%0.3f", (pre[1] + pre[4]) / sum(pre))
  )

  precision <- as.numeric(
    sprintf("%0.3f", pre[4] / (pre[4] + pre[3]))
  )
  recall <- as.numeric(
    sprintf("%0.3f", pre[4] / (pre[4] + pre[2]))
  )
  f1_score <- as.numeric(
    sprintf(
      "%0.3f",
      2 * (precision * recall) / (precision + recall)
    )
  )

  pred_edges <- network_table[predictor_binary == 1, c("regulator", "target")]
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
  intersection_size <- length(overlap_edges)
  union_size <- length(union(pred_edge_ids, true_edge_ids))

  si_score <- as.numeric(
    sprintf("%0.3f", intersection_size)
  )
  ji_score <- as.numeric(
    sprintf(
      "%0.3f",
      ifelse(union_size > 0, intersection_size / union_size, 0)
    )
  )

  if (plot) {
    all_genes <- sort(
      unique(
        c(
          network_table$regulator, network_table$target,
          ground_truth$regulator, ground_truth$target
        )
      )
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
      edge_plot_data <- rbind(
        edge_plot_data,
        overlap_data
      )
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

    auroc_plot <- subset(
      fortify(auc_curves),
      curvetype == "ROC"
    ) |>
      ggplot(aes(x = x, y = y)) +
      geom_line(color = line_color, linewidth = line_width) +
      geom_abline(
        slope = 1, intercept = 0,
        color = line_color,
        linetype = "dotted",
        linewidth = line_width
      ) +
      labs(
        title = paste("AUROC:", aurco),
        x = "False positive rate",
        y = "True positive rate"
      ) +
      theme_minimal() +
      coord_fixed()

    auprc_plot <- subset(
      fortify(auc_curves),
      curvetype == "PRC"
    ) |>
      ggplot(aes(x = x, y = y)) +
      geom_line(color = line_color, linewidth = line_width) +
      labs(
        title = paste("AUPRC:", auprc),
        x = "Recall",
        y = "Precision"
      ) +
      theme_minimal() +
      coord_fixed()

    conf_plot <- data.frame(
      Predicted = rep(c("Negative", "Positive"), each = 2),
      Actual = rep(c("Negative", "Positive"), 2),
      count = pre
    ) |>
      ggplot(aes(x = Actual, y = Predicted)) +
      geom_tile(aes(fill = count), color = "gray50") +
      geom_text(aes(label = count), color = "white") +
      scale_fill_gradient(low = "#4ECDC4", high = "#FF6B6B") +
      theme_minimal() +
      labs(
        title = "Confusion Matrix"
      ) +
      theme(legend.position = "none") +
      coord_fixed()

    network_plot <- ggplot() +
      geom_tile(
        data = edge_plot_data,
        aes(
          y = regulator,
          x = target,
          fill = type
        ),
        color = "white",
        width = 0.9,
        height = 0.9
      ) +
      scale_fill_manual(values = c(
        "Predicted Only" = "#FF6B6B",
        "Ground Truth Only" = "#4ECDC4",
        "Overlapping" = "#FFB347"
      )) +
      theme_minimal() +
      labs(
        title = "Network Edge Comparison",
        y = "Regulator",
        x = "Target"
      ) +
      theme(
        axis.text.x = element_text(
          angle = 45,
          hjust = 1,
          vjust = 1,
          size = 7
        ),
        axis.text.y = element_text(size = 7),
        legend.position = "bottom",
        legend.title = element_blank(),
        panel.grid = element_line(color = "gray95"),
        plot.margin = margin(b = 20)
      ) +
      coord_fixed()

    edge_details <- data.frame(
      Category = c(
        "Overlapping (SI)",
        "Predicted Only",
        "Ground Truth Only",
        "Total Predicted",
        "Total Ground Truth"
      ),
      count = c(
        intersection_size,
        length(setdiff(pred_edge_ids, true_edge_ids)),
        length(setdiff(true_edge_ids, pred_edge_ids)),
        length(pred_edge_ids),
        length(true_edge_ids)
      ),
      Type = c(
        "Overlap",
        "Prediction",
        "Ground Truth",
        "Total",
        "Total"
      )
    )

    edge_stats_plot <- ggplot(
      edge_details,
      aes(
        x = factor(Category,
          levels = c(
            "Total Predicted",
            "Total Ground Truth",
            "Overlapping (SI)",
            "Predicted Only",
            "Ground Truth Only"
          )
        ),
        y = count,
        fill = Type
      )
    ) +
      geom_bar(stat = "identity", width = 0.7) +
      geom_text(aes(label = count),
        vjust = -0.3
      ) +
      scale_fill_manual(values = c(
        "Overlap" = "#FFB347",
        "Prediction" = "#FF6B6B",
        "Ground Truth" = "#4ECDC4",
        "Total" = "#6C757D"
      )) +
      theme_minimal() +
      labs(
        title = "Edge Distribution",
        x = "",
        y = "Count"
      ) +
      theme(
        axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = "bottom",
        legend.title = element_blank()
      ) +
      scale_y_continuous(
        expand = expansion(mult = c(0, 0.1)),
        limits = function(x) c(0, max(x) * 1.05)
      )

    ratio_metrics <- data.frame(
      Metric = c("AUROC", "AUPRC", "ACC", "F1", "JI", "Precision", "Recall"),
      Value = c(
        aurco, auprc, acc, f1_score, ji_score, precision, recall
      )
    )

    metrics_plot <- ggplot(
      ratio_metrics,
      aes(x = stats::reorder(Metric, -Value), y = Value)
    ) +
      geom_bar(stat = "identity", fill = line_color, alpha = 0.8) +
      geom_text(aes(label = Value),
        vjust = -0.5
      ) +
      theme_minimal() +
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

    final_plot <- patchwork::wrap_plots(
      list(
        auroc_plot, auprc_plot, metrics_plot,
        conf_plot, network_plot, edge_stats_plot
      ),
      ncol = 3,
      nrow = 2,
      widths = c(1, 1, 1),
      heights = c(1, 1)
    )

    print(final_plot)
  }

  if (type == "all") {
    return(
      data.frame(
        AUROC = aurco,
        AUPRC = auprc,
        ACC = acc,
        Precision = precision,
        Recall = recall,
        F1 = f1_score,
        JI = ji_score,
        SI = si_score
      )
    )
  } else if (type == "auroc") {
    return(paste("AUROC:", aurco))
  } else if (type == "auprc") {
    return(paste("AUPRC:", auprc))
  } else if (type == "acc") {
    return(paste("ACC:", acc))
  } else if (type == "precision") {
    return(paste("Precision:", precision))
  } else if (type == "recall") {
    return(paste("Recall:", recall))
  } else if (type == "f1") {
    return(paste("F1:", f1_score))
  } else if (type == "si") {
    return(paste("SI:", si_score))
  } else if (type == "ji") {
    return(paste("JI:", ji_score))
  } else {
    stop("Invalid metric type")
  }
}
