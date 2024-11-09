#' @title Calculate metrics
#'
#' @param network_table The weight data table of network
#' @param ground_truth Ground truth for calculate AUC
#' @param type The type of metrics to calculate
#' @param plot If true, draw and print figure of AUC
#' @param line_color The color of line in the figure
#' @param line_width The width of line in the figure
#'
#' @return AUC values and figure
#' @export
#'
#' @examples
#' data("example_matrix")
#' data("example_ground_truth")
#' inferCSN(example_matrix) |>
#'   calculate_metrics(example_ground_truth, plot = TRUE)
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
  acc <- sprintf("%0.3f", (pre[1] + pre[4]) / sum(pre))

  precision <- sprintf("%0.3f", pre[4] / (pre[4] + pre[3]))
  recall <- sprintf("%0.3f", pre[4] / (pre[4] + pre[2]))
  f1_score <- sprintf(
    "%0.3f",
    2 * (as.numeric(precision) * as.numeric(recall)) / (as.numeric(precision) + as.numeric(recall))
  )

  auc_metric <- data.frame(
    AUROC = sprintf("%0.3f", auc$aucs[1]),
    AUPRC = sprintf("%0.3f", auc$aucs[2]),
    ACC = acc,
    Precision = precision,
    Recall = recall,
    F1 = f1_score
  )

  if (plot) {
    plots <- list()

    if (type %in% c("all", "auroc")) {
      auroc <- subset(
        fortify(auc_curves),
        curvetype == "ROC"
      ) |>
        ggplot(
          aes(x = x, y = y)
        ) +
        geom_line(
          color = line_color,
          linewidth = line_width
        ) +
        geom_abline(
          slope = 1,
          intercept = 0,
          color = line_color,
          linetype = "dotted",
          linewidth = line_width
        ) +
        labs(
          title = paste("AUROC:", auc_metric$AUROC),
          x = "False positive rate",
          y = "True positive rate"
        ) +
        xlim(0, 1) +
        ylim(0, 1) +
        coord_fixed() +
        theme_bw()
      plots <- c(plots, list(auroc))
    }

    if (type %in% c("all", "auprc")) {
      auprc <- subset(
        fortify(auc_curves),
        curvetype == "PRC"
      ) |>
        ggplot(
          aes(x = x, y = y)
        ) +
        geom_line(
          color = line_color,
          linewidth = line_width
        ) +
        labs(
          title = paste("AUPRC:", auc_metric$AUPRC),
          x = "Recall",
          y = "Precision"
        ) +
        xlim(0, 1) +
        ylim(0, 1) +
        coord_fixed() +
        theme_bw()
      plots <- c(plots, list(auprc))
    }

    conf_plot <- data.frame(
      Predicted = rep(c("Negative", "Positive"), each = 2),
      Actual = rep(c("Negative", "Positive"), 2),
      count = pre
    ) |>
      ggplot(
        aes(x = Actual, y = Predicted)
      ) +
      geom_tile(
        aes(fill = count),
        color = "gray50"
      ) +
      geom_text(
        aes(label = count),
        color = "gray50"
      ) +
      scale_fill_gradient(
        low = "#69bbd6", high = line_color
      ) +
      labs(
        x = "Actual",
        y = "Predicted"
      ) +
      theme_minimal() +
      theme(
        legend.position = "none"
      ) +
      coord_fixed()
    plots <- c(plots, list(conf_plot))

    if (length(plots) > 0) {
      p <- patchwork::wrap_plots(
        plots,
        ncol = length(plots)
      )
      print(p)
    }
  }

  if (type == "auroc") {
    return(paste("AUROC:", auc_metric$AUROC))
  } else if (type == "auprc") {
    return(paste("AUPRC:", auc_metric$AUPRC))
  } else if (type == "acc") {
    return(paste("ACC:", acc))
  } else if (type == "precision") {
    return(paste("Precision:", precision))
  } else if (type == "recall") {
    return(paste("Recall:", recall))
  } else if (type == "f1") {
    return(paste("F1:", f1_score))
  } else {
    return(auc_metric)
  }
}
