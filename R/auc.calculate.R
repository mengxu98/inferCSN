#' @title AUC value calculate
#'
#' @param weightDT The weight data table of network
#' @param groundTruth Ground truth for calculate AUC
#' @param plot If true, draw and print figure of AUC
#' @param lineColor The color of line in the figure
#' @param lineWidth The width of line in the figure
#'
#' @import patchwork
#' @import ggplot2
#'
#' @return AUC values and figure
#' @export
#'
#' @examples
#' library(inferCSN)
#' data("exampleMatrix")
#' data("exampleGroundTruth")
#' weightDT <- inferCSN(exampleMatrix)
#' auc <- auc.calculate(weightDT, exampleGroundTruth, plot = TRUE)
#' head(auc)
#'
auc.calculate <- function(weightDT,
                          groundTruth,
                          plot = FALSE,
                          lineColor = "#1563cc",
                          lineWidth = 1) {
  # Check input data
  colnames(weightDT) <- c("regulator", "target", "weight")
  weightDT$weight <- abs(as.numeric(weightDT$weight))

  if (ncol(groundTruth) > 2) groundTruth <- groundTruth[, 1:2]
  names(groundTruth) <- c("regulator", "target")
  groundTruth$label <- rep(1, nrow(groundTruth))

  gold <- merge(weightDT, groundTruth,
                by = c("regulator", "target"),
                all.x = TRUE)
  gold$label[is.na(gold$label)] <- 0

  aucCurves <- precrec::evalmod(scores = gold$weight,
                                labels = gold$label)

  auc <- attr(aucCurves, "auc")

  aucMetric <- data.frame(AUROC = rep(0.000, 1),
                          AUPRC = rep(0.000, 1),
                          ACC = rep(0.000, 1))
  aucMetric[1, "AUROC"] <- sprintf("%0.3f", auc$aucs[1])
  aucMetric[1, "AUPRC"] <- sprintf("%0.3f", auc$aucs[2])
  aucMetric[1, "ACC"] <- sprintf("%0.3f", acc.calculate(gold))
  if (plot) {
    # Separate data
    aurocDf <- subset(fortify(aucCurves),
                      curvetype == "ROC")
    auprcDf <- subset(fortify(aucCurves),
                      curvetype == "PRC")

    # Plot
    auroc <- ggplot(aurocDf, aes(x = x, y = y)) +
      geom_line(color = lineColor,
                linewidth = lineWidth) +
      geom_abline(slope = lineWidth,
                  color = lineColor,
                  linetype = "dotted",
                  linewidth = lineWidth) +
      labs(title = paste("AUROC:", aucMetric[1]),
           x = "False positive rate",
           y = "True positive rate") +
      xlim(0, 1) +
      ylim(0, 1) +
      coord_fixed() +
      theme_bw()

    auprc <- ggplot(auprcDf, aes(x = x, y = y)) +
      geom_line(color = lineColor,
                linewidth = 1) +
      labs(title = paste("AUPRC:", aucMetric[2]),
           x = "Recall",
           y = "Precision") +
      xlim(0, 1) +
      ylim(0, 1) +
      coord_fixed() +
      theme_bw()

    # Combine two plots by `patchwork` package
    p <- auroc + auprc
    print(p)
  }

  return(aucMetric)
}

#' ACC calculate
#'
#' @param gold Data
#'
#' @return ACC
#' @export
#'
acc.calculate <- function(gold) {
  results <- pROC::roc(gold$label ~ gold$weight,
                       direction = "<",
                       levels = c(0, 1))

  # After this operation, '0' indicate positive
  reverseLabel <- 2 - as.numeric(as.factor(gold$label))

  sensitivities <- results$sensitivities
  specificities <- results$specificities
  selectValue <- sensitivities + specificities - 1
  cutValueResults <- results$thresholds[selectValue == max(selectValue)]
  selectSensitivities <- sensitivities[selectValue == max(selectValue)]

  cutValue <- cutValueResults[selectSensitivities == max(selectSensitivities)]

  predictorBinary <- rep(0, length(results$predictor))
  predictorBinary[results$predictor >= cutValue] <- 1
  predictorBinary <- as.factor(predictorBinary)
  levels(predictorBinary) <- c("0", "1")
  predictorBinary <- factor(predictorBinary, levels = c(1, 0))

  pre <- as.vector(table(predictorBinary, reverseLabel))

  acc <- (pre[1] + pre[4]) / sum(pre)

  return(acc)
}
