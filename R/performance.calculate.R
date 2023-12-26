#' @title AUC value calculate
#'
#' @param weight_table The weight data table of network
#' @param ground_truth Ground truth for calculate AUC
#' @param plot If true, draw and print figure of AUC
#' @param line_color The color of line in the figure
#' @param line_width The width of line in the figure
#'
#' @import patchwork
#' @import ggplot2
#'
#' @return AUC values and figure
#' @export
#'
#' @examples
#' library(inferCSN)
#' data("example_matrix")
#' data("example_ground_truth")
#' weight_table <- inferCSN(example_matrix)
#' auc <- auc.calculate(weight_table, example_ground_truth, plot = TRUE)
#' auc
auc.calculate <- function(
    weight_table,
    ground_truth,
    plot = FALSE,
    line_color = "#1563cc",
    line_width = 1) {

  gold <- prepare.performance.data(
    weight_table,
    ground_truth)

  auc_curves <- precrec::evalmod(
    scores = gold$weight,
    labels = gold$label
  )

  auc <- attr(auc_curves, "auc")

  auc_metric <- data.frame(
    AUROC = rep(0.000, 1),
    AUPRC = rep(0.000, 1)
  )
  auc_metric[1, "AUROC"] <- sprintf("%0.3f", auc$aucs[1])
  auc_metric[1, "AUPRC"] <- sprintf("%0.3f", auc$aucs[2])
  if (plot) {
    # Separate data
    auroc_table <- subset(
      fortify(auc_curves),
      curvetype == "ROC"
    )
    auprc_table <- subset(
      fortify(auc_curves),
      curvetype == "PRC"
    )

    # Plot
    auroc <- ggplot(auroc_table, aes(x = x, y = y)) +
      geom_line(
        color = line_color,
        linewidth = line_width
      ) +
      geom_abline(
        slope = line_width,
        color = line_color,
        linetype = "dotted",
        linewidth = line_width
      ) +
      labs(
        title = paste("AUROC:", auc_metric[1]),
        x = "False positive rate",
        y = "True positive rate"
      ) +
      xlim(0, 1) +
      ylim(0, 1) +
      coord_fixed() +
      theme_bw()

    auprc <- ggplot(auprc_table, aes(x = x, y = y)) +
      geom_line(
        color = line_color,
        linewidth = 1
      ) +
      labs(
        title = paste("AUPRC:", auc_metric[2]),
        x = "Recall",
        y = "Precision"
      ) +
      xlim(0, 1) +
      ylim(0, 1) +
      coord_fixed() +
      theme_bw()

    # Combine two plots by `patchwork` package
    p <- auroc + auprc
    print(p)
  }

  return(auc_metric)
}

#' ACC calculate
#'
#' @inheritParams auc.calculate
#'
#' @return ACC value
#' @export
#'
#' @examples
#' library(inferCSN)
#' data("example_matrix")
#' data("example_ground_truth")
#' weight_table <- inferCSN(example_matrix)
#' acc <- acc.calculate(weight_table, example_ground_truth)
#' acc
acc.calculate <- function(
    weight_table,
    ground_truth) {

  gold <- prepare.performance.data(
    weight_table,
    ground_truth)
  results <- pROC::roc(
    gold$label ~ gold$weight,
    direction = "<",
    levels = c(0, 1))

  # After this operation, '0' indicate positive
  reverse_label <- 2 - as.numeric(as.factor(gold$label))

  sensitivities <- results$sensitivities
  specificities <- results$specificities
  select_value <- sensitivities + specificities - 1
  cut_value_results <- results$thresholds[select_value == max(select_value)]
  select_sensitivities <- sensitivities[select_value == max(select_value)]

  cut_value <- cut_value_results[
    select_sensitivities == max(select_sensitivities)
  ]

  predictor_binary <- rep(0, length(results$predictor))
  predictor_binary[results$predictor >= cut_value] <- 1
  predictor_binary <- as.factor(predictor_binary)
  levels(predictor_binary) <- c("0", "1")
  predictor_binary <- factor(predictor_binary, levels = c(1, 0))

  pre <- as.vector(table(predictor_binary, reverse_label))

  acc <- (pre[1] + pre[4]) / sum(pre)
  acc <- sprintf("%0.3f", acc)
  return(acc)
}

#' @title prepare.performance.data
#'
#' @inheritParams auc.calculate
#'
#' @return Formated data
#' @export
prepare.performance.data <- function(
    weight_table,
    ground_truth) {
  # Check input data
  colnames(weight_table) <- c("regulator", "target", "weight")
  weight_table$weight <- abs(as.numeric(weight_table$weight))

  if (ncol(ground_truth) > 2) ground_truth <- ground_truth[, 1:2]
  names(ground_truth) <- c("regulator", "target")
  ground_truth$label <- rep(1, nrow(ground_truth))

  gold <- merge(
    weight_table,
    ground_truth,
    by = c("regulator", "target"),
    all.x = TRUE
  )
  gold$label[is.na(gold$label)] <- 0
  return(gold)
}
