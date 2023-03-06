#' inferCSN
#'
#' @param matrix Gene experssion matrix
#' @param penalty Default "L0"
#' @param regulators Regulator genes
#' @param targets Target genes
#' @param maxSuppSize The number of non-zore coef
#' @param cores CPU cores
#'
#' @return A list of gene-gene regulatory relationship
#' @export
#'
#' @examples
#'   data("exampleDataMatrix")
#'   weightList <- inferCSN(exampleDataMatrix)
inferCSN <- function(matrix,
                    penalty = NULL,
                    regulators = NULL,
                    targets = NULL,
                    maxSuppSize = NULL,
                    cores = 1) {
  if (!require("L0Learn")) devtools::install_github("hazimehh/L0Learn")
  matrix <- as.data.frame(matrix)

  if (is.null(penalty)) penalty <- "L0"

  if (is.null(maxSuppSize)) maxSuppSize <- dim(matrix)[2]

  if (is.null(targets)) targets <- colnames(matrix)

  if (!is.null(regulators)) {
    matrix <- matrix[, regulators]
  } else {
    regulators <- colnames(matrix)
  }

  if (cores == 1) {
    weightList <- c()
    for (i in 1:length(regulators)) {
      message(paste("Running for", i, "of", length(regulators), "gene:", regulators[i]))
      X <- as.matrix(matrix[, -which(colnames(matrix) == regulators[i])])
      Y <- matrix[, regulators[i]]
      temp <- inferCSN.fit(X, Y,
                     penalty = penalty,
                     nFolds = 10,
                     seed = 1,
                     maxSuppSize = maxSuppSize,
                     nGamma = 5,
                     gammaMin = 0.0001,
                     gammaMax = 10
      )
      temp <- as.vector(temp)
      wghts <- temp[-1]
      wghts <- abs(wghts)
      wghts <- wghts / sum(wghts)
      if (F) {
        wghts <- wghts / max(wghts)
        indices <- sort.list(wghts, decreasing = TRUE)
        zeros <- which(wghts <= 0.8)
        # wghts[1:length(wghts)] <- 1
        wghts[zeros] <- 0
      }
      if (length(wghts) != dim(X)[2]) {
        weightd <- data.frame(regulatoryGene = colnames(X), targetGene = regulators[i], weight = 0)
      } else {
        weightd <- data.frame(regulatoryGene = colnames(X), targetGene = regulators[i], weight = wghts)
      }
      weightList <- rbind.data.frame(weightList, weightd)
      if (i == length(regulators)) {
        weightList <- weightList[order(weightList$weight, decreasing = TRUE), ]
      }
    }
  } else {
    cores <- min(parallel::detectCores(logical = F), cores)
    cl <- parallel::makeCluster(cores)
    doParallel::registerDoParallel(cl)
    # doParallel::registerDoParallel(cores = cores)
    message(paste("\nUsing", foreach::getDoParWorkers(), "cores."))
    "%dopar%" <- foreach::"%dopar%"
    suppressPackageStartupMessages(
      weightList <- doRNG::"%dorng%"(foreach::foreach(regulator = regulators, .combine = "rbind", .export = "inferCSN.fit"), {
        X <- as.matrix(matrix[, -which(colnames(matrix) == regulator)])
        Y <- matrix[, regulator]
        temp <- inferCSN.fit(X, Y,
                       penalty = penalty,
                       nFolds = 10,
                       seed = 1,
                       maxSuppSize = maxSuppSize,
                       nGamma = 5,
                       gammaMin = 0.0001,
                       gammaMax = 10
        )
        temp <- as.vector(temp)
        wghts <- temp[-1]
        wghts <- abs(wghts)
        wghts <- wghts / sum(wghts)

        if (length(wghts) != dim(X)[2]) {
          weightd <- data.frame(regulatoryGene = colnames(X), targetGene = regulator, weight = 0)
        } else {
          weightd <- data.frame(regulatoryGene = colnames(X), targetGene = regulator, weight = wghts)
        }
      })
    )
    parallel::stopCluster(cl)
  }
  attr(weightList, "rng") <- NULL
  attr(weightList,"doRNG_version") <- NULL
  weightList <- weightList[order(weightList$weight, decreasing = TRUE), ]
  return(weightList)
}
