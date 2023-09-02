#' @title Inferring cell-specific gene regulatory network
#'
#' @param matrix An expression matrix, cells by genes
#' @param penalty [Default = "L0"] The type of regularization.
#' This can take either one of the following choices: "L0"and "L0L2".
#' For high-dimensional and sparse data, such as single-cell transcriptome data, "L0L2" is more effective
#' @param algorithm [Default = "CD"] Currently "CD" and "CDPSI" are supported.
#' The CDPSI algorithm may yield better results, but it also increases running time
#' @param crossValidation [Default = FALSE] Check whether cross validation is used
#' @param nFolds [Default = 10] The number of folds for cross-validation
#' @param regulators [Default = NULL] Regulator genes
#' @param targets [Default = NULL] Target genes
#' @param maxSuppSize [Default = NULL] The number of non-zore coef, this value will affect the final performance
#' @param verbose [Default = FALSE] Print detailed information
#' @param cores [Default = 1] CPU cores
#'
#' @import magrittr
#' @importFrom utils methods read.table setTxtProgressBar txtProgressBar
#'
#' @return A data table of gene-gene regulatory relationship
#' @export
#'
#' @examples
#' library(inferCSN)
#' data("exampleMatrix")
#' weightDT <- inferCSN(exampleMatrix, verbose = TRUE)
#' head(weightDT)
#' weightDT <- inferCSN(exampleMatrix, cores = 2)
#' head(weightDT)
#'
inferCSN <- function(matrix,
                     penalty = "L0",
                     algorithm = "CD",
                     crossValidation = FALSE,
                     nFolds = 10,
                     regulators = NULL,
                     targets = NULL,
                     maxSuppSize = NULL,
                     verbose = FALSE,
                     cores = 1) {
  # Check the penalty term of the regression model
  if (!any(c("L0", "L0L2") == penalty)) {
    stop("inferCSN does not support '", penalty, "' penalty regression......\n",
         "Please set penalty item as 'L0' or 'L0L2'......")
  }

  # Check the algorithm of the regression model
  if (!any(c("CD", "CDPSI") == algorithm)) {
    stop("inferCSN does not support '", algorithm, "' algorithm......\n",
         "Please set algorithm as 'CD' or 'CDPSI'......")
  }

  if (!is.null(regulators)) {
    regulatorsMatrix <- matrix[, intersect(colnames(matrix), regulators)]
  } else {
    regulatorsMatrix <- matrix
  }
  regulators <- colnames(regulatorsMatrix)

  if (!is.null(targets)) {
    targetsMatrix <- matrix[, intersect(colnames(matrix), targets)]
  } else {
    targetsMatrix <- matrix
  }
  targets <- colnames(targetsMatrix)

  cores <- min(parallel::detectCores(logical = FALSE), cores, length(targets))
  if (cores == 1) {
    # Format progress information
    pb <- progress::progress_bar$new(format = "Running [:bar] :percent, No.:current of :total gene,:elapsed......",
                                     total = length(targets),
                                     clear = TRUE,
                                     width = 100)

    weightDT <- c()
    for (i in 1:length(targets)) {
      if (verbose) {
        # Print progress
        pb$tick()
        Sys.sleep(0.05)
      }

      weightDT <- rbind(weightDT,
                        sub.inferCSN(regulatorsMatrix = regulatorsMatrix,
                                     targetsMatrix = targetsMatrix,
                                     target = targets[i],
                                     crossValidation = crossValidation,
                                     penalty = penalty,
                                     algorithm = algorithm,
                                     nFolds = nFolds,
                                     maxSuppSize = maxSuppSize,
                                     verbose = verbose))
    }

  } else {
    cl <- snow::makeSOCKcluster(cores)
    doSNOW::registerDoSNOW(cl)
    pb <- txtProgressBar(min = 1,
                         max = length(targets),
                         width = 100,
                         style = 3)

    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    "%dopar%" <- foreach::"%dopar%"
    if (verbose) {
      multipleCoresSetting <- foreach::foreach(target = targets,
                                               .combine = "rbind",
                                               .export = c("sub.inferCSN", "sparse.regression"),
                                               .options.snow = opts)
    } else {
      multipleCoresSetting <- foreach::foreach(target = targets,
                                               .combine = "rbind",
                                               .export = c("sub.inferCSN", "sparse.regression"))
    }

    weightDT <- multipleCoresSetting %dopar% {
      sub.inferCSN(regulatorsMatrix = regulatorsMatrix,
                   targetsMatrix = targetsMatrix,
                   target = target,
                   crossValidation = crossValidation,
                   penalty = penalty,
                   algorithm = algorithm,
                   nFolds = nFolds,
                   maxSuppSize = maxSuppSize,
                   verbose = verbose)
    }

    if (verbose) close(pb)
    parallel::stopCluster(cl)
  }

  weightDT <- weightDT[order(abs(weightDT$weight), decreasing = TRUE), ]
  return(weightDT)
}
