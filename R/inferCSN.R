#' inferCSN
#'    A method for inferring cell-specific gene regulatory network from single-cell transcriptome data.
#'
#' @param data An object
#' @param normalize Data normalize
#' @param penalty Default "L0"
#' @param crossValidation Cross validation
#' @param nFolds N folds cross validation
#' @param regulators Regulator genes
#' @param targets Target genes
#' @param maxSuppSize The number of non-zore coef
#' @param verbose Print detailed information
#' @param cores CPU cores
#' @param nGamma nGamma
#'
#' @return A list of gene-gene regulatory relationship
#' @export
#'
#' @examples
#'   data("exampleDataMatrix")
#'   weightList <- inferCSN(exampleDataMatrix)
inferCSN <- function(data = NULL,
                     normalize = FALSE,
                     penalty = NULL,
                     crossValidation = FALSE,
                     nFolds = 10,
                     regulators = NULL,
                     targets = NULL,
                     maxSuppSize = NULL,
                     nGamma = 5,
                     verbose = FALSE,
                     cores = 1) {
  # Processing
  if (!is.null(data)) {
    matrix <- data.processing(data, normalize = normalize, verbose = verbose)
  } else {
    stop("Please ensure provide an object!")
  }

  if (is.null(penalty)) penalty <- "L0"

  if (is.null(maxSuppSize)) maxSuppSize <- ncol(matrix)

  # if (is.null(targets)) targets <- colnames(matrix)

  if (!is.null(targets)) {
    targetsMatrix <- matrix[, targets]
  } else {
    targets <- colnames(matrix)
    targetsMatrix <- matrix
  }

  # if (!is.null(regulators)) regulators <- colnames(matrix)

  if (!is.null(regulators)) {
    regulatorsMatrix <- matrix[, regulators]
  } else {
    regulators <- colnames(matrix)
    regulatorsMatrix <- matrix
  }

  if (cores == 1) {
    weightList <- c()
    for (i in 1:length(targets)) {
      if (verbose) message(paste("Running for", i, "of", length(targets), "gene:", targets[i],"......"))

      if (targets[i] %in% regulators) {
        X <- as.matrix(targetsMatrix[, -which(colnames(targetsMatrix) == targets[i])])
      } else {
        X <- as.matrix(targetsMatrix)
      }

      y <- regulatorsMatrix[, targets[i]]

      temp <- inferCSN.fit(X, y,
                           penalty = penalty,
                           crossValidation = crossValidation,
                           nFolds = nFolds,
                           maxSuppSize = maxSuppSize,
                           nGamma = nGamma,
                           verbose = verbose
      )
      temp <- as.vector(temp)
      wghts <- temp[-1]
      wghts <- abs(wghts)
      wghts <- wghts / sum(wghts)
      if (length(wghts) != ncol(X)) {
        weightd <- data.frame(regulator = colnames(X), target = targets[i], weight = 0)
      } else {
        weightd <- data.frame(regulator = colnames(X), target = targets[i], weight = wghts)
      }
      # if (length(wghts) != ncol(X)) weight <- matrix(0, nrow = ncol(X), ncol = 1)
      # weightd <- data.frame(regulator = colnames(X), target = targets[i], weight = wghts)
      weightList <- rbind.data.frame(weightList, weightd)
    }
  } else {
    cores <- min(parallel::detectCores(logical = F), cores)
    cl <- parallel::makeCluster(cores)
    doParallel::registerDoParallel(cl)
    # doParallel::registerDoParallel(cores = cores)
    if (verbose) message(paste("\nUsing", foreach::getDoParWorkers(), "cores......"))
    "%dopar%" <- foreach::"%dopar%"
    suppressPackageStartupMessages(
      weightList <- doRNG::"%dorng%"(
        foreach::foreach(target = targets,
                         .combine = "rbind",
                         .export = "inferCSN.fit",
                         .packages = c("doParallel")), {
                           # X <- as.matrix(matrix[, -which(colnames(matrix) == target)])
                           # y <- matrix[, target]

                           if (target %in% regulators) {
                             X <- as.matrix(targetsMatrix[, -which(colnames(targetsMatrix) == target)])
                           } else {
                             X <- as.matrix(targetsMatrix)
                           }

                           y <- regulatorsMatrix[, target]

                           temp <- inferCSN.fit(X, y,
                                                penalty = penalty,
                                                crossValidation = crossValidation,
                                                nFolds = nFolds,
                                                maxSuppSize = maxSuppSize,
                                                nGamma = nGamma,
                                                verbose = verbose
                           )
                           temp <- as.vector(temp)
                           wghts <- temp[-1]
                           wghts <- abs(wghts)
                           wghts <- wghts / sum(wghts)
                           if (length(wghts) != ncol(X)) {
                             weightd <- data.frame(regulator = colnames(X), target = target, weight = 0)
                           } else {
                             weightd <- data.frame(regulator = colnames(X), target = target, weight = wghts)
                           }
                         }
        )
    )
    parallel::stopCluster(cl)
  }
  attr(weightList, "rng") <- NULL
  attr(weightList, "doRNG_version") <- NULL
  weightList <- weightList[order(weightList$weight, decreasing = TRUE), ]
  return(weightList)
}
