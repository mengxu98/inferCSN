#' @title auc.calculate
#' @description AUC value calculate
#'
#' @param weightList weightList
#' @param groundTruth groundTruth
#' @param groundTruthTPath groundTruthTPath
#' @param plot [Default = FALSE] plot
#'
#' @import patchwork
#'
#' @return AUC values and AUC plot
#' @export
#'
#' @examples
#' data("exampleDataMatrix")
#' data("exampleDataGroundTruth")
#' weightList <- inferCSN(exampleDataMatrix, cores = 2)
#' auc <- auc.calculate(weightList, exampleDataGroundTruth)
auc.calculate <- function(weightList = NULL,
                          groundTruth = NULL,
                          groundTruthTPath = NULL,
                          plot = FALSE) {
  # Check input data
  if (!is.null(weightList) & ncol(weightList) == 3) {
    names(weightList) <- c("regulator", "target", "weight")
  } else {
    stop("Please provide a data table of regulatory relationships......")
  }

  if (!is.null(groundTruth)) {
    gold_ref <- groundTruth
  } else if (!is.null(groundTruthTPath)) {
    gold_ref <- read.table(
      file = groundTruthTPath,
      header = TRUE,
      sep = ",",
      stringsAsFactors = FALSE
    )
  } else {
    stop("Please provide a reference data or path of ground-truth......")
  }

  if (ncol(gold_ref) > 2) gold_ref <- gold_ref[, 1:2]
  names(gold_ref) <- c("regulator", "target")

  auc.metric <- data.frame(AUROC = rep(0.000, 1), AUPRC = rep(0.000, 1))
  gold_ref$gold <- rep(1, nrow(gold_ref))
  gold <- merge(weightList, gold_ref, by = c("regulator", "target"), all.x = TRUE)
  gold$gold[is.na(gold$gold)] <- 0
  auc.curves <- precrec::evalmod(scores = gold$weight, labels = gold$gold)

  auc <- attr(auc.curves, "auc")
  auc.metric[1, "AUROC"] <- auc$aucs[1]
  auc.metric[1, "AUPRC"] <- auc$aucs[2]

  if (plot) {
    # p <- ggplot2::autoplot(auc.curves)

    # Subset data to separate prc and roc
    auprcDf <- subset(ggplot2::fortify(auc.curves), curvetype == "PRC")
    aurocDf <- subset(ggplot2::fortify(auc.curves), curvetype == "ROC")

    # AUROC
    auroc <- ggplot2::ggplot(aurocDf, ggplot2::aes(x = x, y = y)) +
      ggplot2::geom_line() +
      ggplot2::geom_abline(slope = 1, color = "gray", linetype = "dotted") +
      ggplot2::ggtitle(paste("AUROC:", auc$aucs[1])) +
      ggplot2::labs(x = "FPR", y = "TPR") +
      ggplot2::coord_fixed() +
      ggplot2::theme_bw()

    # AUPRC
    auprc <- ggplot2::ggplot(auprcDf, ggplot2::aes(x = x, y = y)) +
      ggplot2::geom_line() +
      ggplot2::geom_hline(yintercept = 0.5, color = "gray", linetype = "dotted") +
      ggplot2::ggtitle(paste("AUPRC:", auc$aucs[2])) +
      ggplot2::labs(x = "TPR", y = "PPV") +
      ggplot2::ylim(0, 1) +
      ggplot2::coord_fixed() +
      ggplot2::theme_bw()

    # Combine two plots by patchwork
    p <- auroc + auprc

    print(p)
    # Save figure
    cowplot::ggsave2(file = "AUC.pdf", p, width = 18, height = 10, units = "cm")
  }
  return(auc.metric)
}
