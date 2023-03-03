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
#'   weightList <- inferCSN(exampleDataMatrix)
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

  if (ncol(gold_ref) > 2) {
    gold_ref <- gold_ref[, 1:2]
    names(gold_ref) <- c("regulator", "target")
  }else{
    names(gold_ref) <- c("regulator", "target")
  }

  gold_ref$gold <- rep(1, dim(gold_ref)[1])
  gold <- merge(weightList, gold_ref, by = c("regulator", "target"), all.x = TRUE)
  gold$gold[is.na(gold$gold)] <- 0

  auc.metric <- data.frame(AUROC = rep(0.0, 1), AUPRC = rep(0.0, 1))

  auc.metric[1, "AUROC"] <- auroc.curve(obs = gold$gold, pred = gold$weight, plot = plot)
  auc.metric[1, "AUPRC"] <- auprc.curve(obs = gold$gold, pred = gold$weight, plot = plot)

  return(auc.metric)
}

#' auroc.curve
#'  auroc
#' @param pred pred
#' @param obs obs
#' @param plot plot
#'
#' @return auroc
#' @export
#'
auroc.curve <- function(pred, obs, plot = plot) {
  predict <- ROCR::prediction(pred, obs)
  curve <- ROCR::performance(predict, "tpr", "fpr")
  auroc <- ROCR::performance(predict, "auc")@y.values[[1]]

  auroc_data <- data.frame(
    tpr_min = curve@y.values[[1]],
    fpr_min = curve@x.values[[1]],
    auroc = auroc
  )

  if (plot) {
    requireNamespace("ggplot2")
    p <- ggplot() +
      geom_line(data = auroc_data, aes(x = fpr_min, y = tpr_min), color = "#003399", size = 1) +
      geom_line(aes(x = c(0, 1), y = c(0, 1)), color = "grey", size = 1, linetype = 6) +
      theme_bw() +
      annotate("text",
               x = .2, y = .9,
               label = paste("AUROC =", round(auroc_data$auroc, 3)), color = "#003399"
      ) +
      scale_x_continuous(name = "False Positive Rate") +
      scale_y_continuous(name = "True Positive Rate") +
      ggtitle("AUROC")
    print(p)
  }

  return(auroc)
}

#' auprc.curve
#'  auprc
#' @param pred pred
#' @param obs obs
#' @param plot plot
#'
#' @return auprc
#' @export
#'
auprc.curve <- function(pred, obs, plot = plot) {
  predict <- modEvA::AUC(obs = obs, pred = pred, curve = "PR", interval = 0.001, plot = FALSE) #, simplif = TRUE
  # predict <- AUC(obs = obs, pred = pred, trapezoid = "PR", interval = 0.001, plot = FALSE) #, simplif = TRUE
  auprc <- predict$AUC
  auprc_data <- data.frame(
    precision = predict$thresholds$precision,
    recall = predict$thresholds$sensitivity,
    auprc = auprc
  )

  if (plot) {
    requireNamespace("ggplot2")
    p <- ggplot() +
      geom_line(data = auprc_data, aes(x = recall, y = precision), color = "#f97e11", size = 1) +
      geom_line(aes(x = c(0, 1), y = c(0, 1)), color = "grey", size = 1, linetype = 6) +
      theme_bw() +
      annotate("text",
               x = .2, y = .9,
               label = paste("AUPRC =", round(auprc_data$auprc, 3)), color = "#f97e11"
      ) +
      scale_x_continuous(name = "Recall") +
      scale_y_continuous(name = "Precision") +
      ggtitle("Precision-Recall Curve")
    print(p)
  }

  return(auprc)
}
