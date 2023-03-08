#' auc.calculate
#'  AUC value calculate
#'
#' @param weightList weightList
#' @param goldRef goldRef
#' @param goldPathway goldPathway
#' @param plot plot
#'
#' @return AUC
#' @export
#'
#' @examples
#'   data("exampleDataMatrix")
#'   data("exampleDataGroundTruth")
#'   weightList <- inferCSN(exampleDataMatrix, cores = 2)
#'   auc <- auc.calculate(weightList, exampleDataGroundTruth)
auc.calculate <- function(weightList = NULL,
                          goldRef = NULL,
                          goldPathway = NULL,
                          plot = FALSE) {
  if (!is.null(weightList) & ncol(weightList) == 3) {
    names(weightList) <- c("regulator", "target", "weight")
  } else {
    stop("Please provide a 'weightList' of gene-gene relationships!")
  }

  if (!is.null(goldRef) | !is.null(goldPathway)) {
    if (!is.null(goldRef)) {
      gold_ref <- goldRef
    } else {
      gold_ref <- read.table(
        file = goldPathway,
        stringsAsFactors = FALSE,
        header = TRUE,
        sep = ","
      )
    }
  } else {
    stop("Please provide a reference dataset or goldPathway of gene-gene relationships!")
  }

  if (ncol(gold_ref) > 2) gold_ref <- gold_ref[, 1:2]
  names(gold_ref) <- c("regulator", "target")

  auc.metric <- data.frame(AUROC = rep(0.000, 1), AUPRC = rep(0.000, 1))
  gold_ref$gold <- rep(1, dim(gold_ref)[1])
  gold <- merge(weightList, gold_ref, by = c("regulator", "target"), all.x = TRUE)
  gold$gold[is.na(gold$gold)] <- 0
  auc.curves <- precrec::evalmod(scores = gold$weight, labels = gold$gold)
  if (plot) {
    p <- ggplot2::autoplot(auc.curves)
    cowplot::ggsave2(file = "ROC and Precision-Recall.png", p, width = 18, height = 10, units = "cm")
  }
  auc <- attr(auc.curves, "auc")
  auc.metric[1, "AUROC"] <- auc$aucs[1]
  auc.metric[1, "AUPRC"] <- auc$aucs[2]
  return(auc.metric)
}
