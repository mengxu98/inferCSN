#' Create GRN table
#'
#' @param zscores zscores
#' @param corrMatrix zscores
#' @param genes zscores
#' @param threshold zscores
#' @param self_remove zscores
#'
#' @return data frame of GRN
#'
#' @export
extract_net <- function(
    zscores,
    corrMatrix,
    genes,
    threshold,
    self_remove = TRUE) {
  targets <- vector()
  regulators <- vector()
  zscoresX <- vector()
  correlations <- vector()

  # targets <- rep("", 1e6)
  # regulators <- rep("", 1e6)
  # zscoresX <- rep(0, 1e6)
  # correlations <- rep(0, 1e6)

  str <- 1
  stp <- 1
  for (target in genes) {
    x <- zscores[target, ]
    regs <- names(which(x > threshold))
    if (length(regs) > 0) {
      zzs <- x[regs]
      corrs <- corrMatrix[target, regs]
      ncount <- length(regs)
      stp <- str + ncount - 1
      targets[str:stp] <- rep(target, ncount)
      regulators[str:stp] <- regs
      zscoresX[str:stp] <- zzs
      correlations[str:stp] <- corrs
      str <- stp + 1
    }
  }
  targets <- targets[1:stp]
  regulators <- regulators[1:stp]
  zscoresX <- zscoresX[1:stp]
  correlations <- correlations[1:stp]

  grn <- data.frame(
    target = targets,
    reg = regulators,
    zscore = zscoresX,
    corr = correlations
  )
  colnames(grn)[1:2] <- c("TG", "TF")

  if (self_remove) {
    tfx <- as.vector(grn$TF)
    tgx <- as.vector(grn$TG)
    xi <- which(tfx != tgx)
    grn <- grn[xi, ]
  }
  grn
}
