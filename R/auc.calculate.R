#' @title AUC value calculate
#'
#' @param weightDT The weight data table of network
#' @param groundTruth Ground truth for calculate AUC
#' @param plot If true, draw and print figure of AUC
#' @param fileSave The figure name
#' @param interaction If true, consider the positivity/negativity of interaction
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
                          fileSave = NULL,
                          interaction = FALSE) {
  # Check input data
  colnames(weightDT) <- c("regulator", "target", "weight")
  if (!interaction) weightDT$weight <- abs(as.numeric(weightDT$weight))

  if (ncol(groundTruth) > 2) groundTruth <- groundTruth[, 1:2]
  names(groundTruth) <- c("regulator", "target")

  groundTruth$gold <- rep(1, nrow(groundTruth))
  gold <- merge(weightDT, groundTruth, by = c("regulator", "target"), all.x = TRUE)
  gold$gold[is.na(gold$gold)] <- 0
  aucCurves <- precrec::evalmod(scores = gold$weight, labels = gold$gold)

  auc <- attr(aucCurves, "auc")
  aucMetric <- data.frame(AUROC = rep(0.000, 1), AUPRC = rep(0.000, 1))
  aucMetric[1, "AUROC"] <- sprintf("%0.3f", auc$aucs[1])
  aucMetric[1, "AUPRC"] <- sprintf("%0.3f", auc$aucs[2])

  if (plot) {
    # Subset data to separate prc and roc
    auprcDf <- subset(fortify(aucCurves), curvetype == "PRC")
    aurocDf <- subset(fortify(aucCurves), curvetype == "ROC")

    # Plot
    auroc <- ggplot(aurocDf, aes(x = x, y = y)) +
      geom_line() +
      geom_abline(slope = 1,
                  color = "gray",
                  linetype = "dotted",
                  linewidth = 1) +
      labs(title = paste("AUROC:", aucMetric[1]),
           x = "False positive rate",
           y = "True positive rate") +
      xlim(0, 1) +
      ylim(0, 1) +
      coord_fixed() +
      theme_bw()

    auprc <- ggplot(auprcDf, aes(x = x, y = y)) +
      geom_line() +
      labs(title = paste("AUPRC:", aucMetric[2]),
           x = "Recall",
           y = "Precision") +
      xlim(0, 1) +
      ylim(0, 1) +
      coord_fixed() +
      theme_bw()

    # Combine two plots by patchwork
    p <- auroc + auprc
    print(p)

    # Save figure
    if (!is.null(fileSave)) {
      if (!grepl(".*\\.(pdf|png|jpe?g)$", fileSave)) {
        fileSave <- paste0(fileSave, ".png")
      }
      cowplot::ggsave2(file = fileSave,
                       p,
                       width = 7,
                       height = 3,
                       dpi = 600)
    }
  }
  return(aucMetric)
}
