#' @title auc.calculate
#' @description AUC value calculate
#'
#' @param weightDT The weight data table of network.
#' @param groundTruth Ground truth for calculate AUC.
#' @param plot If true, draw and print figure of AUC.
#' @param fileSave The figure name
#' @param interaction If true, consider the positivity and negativity of interaction.
#'
#' @import patchwork
#'
#' @return AUC values and figure
#' @export
#'
#' @examples
#' data("exampleDataMatrix")
#' data("exampleDataGroundTruth")
#' weightDT <- inferCSN(exampleDataMatrix, cores = 1, verbose = TRUE, algorithm = "CDPSI")
#' auc <- auc.calculate(weightDT, exampleDataGroundTruth, plot = TRUE)
#' head(auc)
auc.calculate <- function(weightDT = NULL,
                          groundTruth = NULL,
                          plot = FALSE,
                          fileSave = NULL,
                          interaction = FALSE) {
  # Check input data
  if (!is.null(weightDT) & ncol(weightDT) == 3) {
    names(weightDT) <- c("regulator", "target", "weight")
    if (!interaction) weightDT$weight <- abs(weightDT$weight)
  } else {
    stop("Please provide the data of regulatory relationships......")
  }

  if (is.null(groundTruth)) stop("Please provide the ground-truth......")
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
    # p <- ggplot2::autoplot(aucCurves)

    # Subset data to separate prc and roc
    auprcDf <- subset(ggplot2::fortify(aucCurves), curvetype == "PRC")
    aurocDf <- subset(ggplot2::fortify(aucCurves), curvetype == "ROC")

    # AUROC plot
    auroc <- ggplot2::ggplot(aurocDf, ggplot2::aes(x = x, y = y)) +
      ggplot2::geom_line() +
      ggplot2::geom_abline(slope = 1, color = "gray", linetype = "dotted") +
      ggplot2::ggtitle(paste("AUROC:", aucMetric[1])) +
      ggplot2::labs(x = "False positive rate", y = "True positive rate") +
      ggplot2::coord_fixed() +
      ggplot2::theme_bw()

    # AUPRC plot
    auprc <- ggplot2::ggplot(auprcDf, ggplot2::aes(x = x, y = y)) +
      ggplot2::geom_line() +
      ggplot2::geom_hline(yintercept = 0.5, color = "gray", linetype = "dotted") +
      ggplot2::ggtitle(paste("AUPRC:", aucMetric[2])) +
      ggplot2::labs(x = "Recall", y = "Precision") +
      ggplot2::ylim(0, 1) +
      ggplot2::coord_fixed() +
      ggplot2::theme_bw()

    # Combine two plots by patchwork
    p <- auroc + auprc
    print(p)

    # Save figure
    if (!is.null(fileSave)) {
      if (figure.format(fileSave)) {
        newFileSave <- fileSave
      } else {
        newFileSave <- paste0(fileSave, ".png")
      }
      cowplot::ggsave2(file = newFileSave,
                       p,
                       width = 18,
                       height = 10,
                       units = "cm",
                       dpi = 600)
    }
  }
  return(aucMetric)
}

#' @title figure.format
#' @description Check figure format
#' @param string A string of file name
#'
#' @return Logic value
#' @export
#'
figure.format <- function(string) {
  logic <- grepl(".*\\.(pdf|png|jpe?g)$", string)
  return(logic)
}
