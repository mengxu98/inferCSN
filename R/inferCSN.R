#' @title inferCSN
#' @description A method for inferring cell-type-specific gene regulatory network
#' from single-cell RNA data.
#'
#' @param matrix An expression matrix, cells by genes.
#' @param penalty [Default = "L0"] The type of regularization.
#' This can take either one of the following choices: "L0"and "L0L2".
#' For high-dimensional and sparse data, such as single-cell transcriptome data, "L0L2" is more effective.
#' @param algorithm [Default = "CD"] Currently "CD" and "CDPSI" are supported.
#' The CDPSI algorithm may yield better results, but it also increases running time.
#' @param crossValidation [Default = FALSE] Check whether cross validation is used.
#' @param nFolds [Default = 10] The number of folds for cross-validation.
#' @param regulators [Default = NULL] Regulator genes.
#' @param targets [Default = NULL] Target genes.
#' @param maxSuppSize [Default = NULL] The number of non-zore coef, this value will affect the final performance.
#' @param verbose [Default = FALSE] Print detailed information.
#' @param cores [Default = 1] CPU cores.
#'
#' @import magrittr
#' @importFrom utils methods read.table setTxtProgressBar txtProgressBar
#'
#' @return A data table of gene-gene regulatory relationship
#' @export
#'
#' @examples
#' data("exampleMatrix")
#' weightDT <- inferCSN(exampleMatrix, verbose = TRUE)
#' head(weightDT)
inferCSN <- function(matrix = NULL,
                     penalty = NULL,
                     algorithm = NULL,
                     crossValidation = FALSE,
                     nFolds = 10,
                     regulators = NULL,
                     targets = NULL,
                     maxSuppSize = NULL,
                     verbose = FALSE,
                     cores = 1) {
  # Data processing
  if (is.null(matrix)) stop("Please ensure provide an expression matrix......")

  # Check the penalty terms of the regression model
  if (!is.null(penalty)) {
    if (!any(c("L0", "L0L2") == penalty)) {
      stop(paste(
        "Note: inferCSN does not support", penalty, "penalty regression......\n",
        "Please set penalty item as 'L0' or 'L0L2'......"
      ))
    }
  } else {
    penalty <- "L0"
  }

  # Check whether cross validation is used
  if (verbose & crossValidation) {
    if (verbose) message(paste("Using", penalty, "penalty and cross validation......"))
  } else {
    if (verbose) message(paste("Using", penalty, "penalty......"))
  }

  # Check the algorithm of the regression model
  if (!is.null(algorithm)) {
    if (!any(c("CD", "CDPSI") == algorithm)) {
      stop(paste(
        "Note: inferCSN does not support", algorithm, "algorithm......\n",
        "Please set algorithm as 'CD' or 'CDPSI'......"
      ))
    }
  } else {
    algorithm <- "CD"
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

  if (cores == 1) {
    if (verbose) {
      # Format progress information
      pb <- progress::progress_bar$new(
        format = "Running [:bar] :percent, No.:current of :total gene,:elapsed......",
        total = length(targets),
        clear = TRUE,
        width = 100
      )
    }

    weightDT <- c()
    for (i in 1:length(targets)) {
      if (verbose) {
        # Print progress
        pb$tick()
        Sys.sleep(0.05)
      }
      subWeightDT <- sub.inferCSN(regulatorsMatrix = regulatorsMatrix,
                                      targetsMatrix = targetsMatrix,
                                      target = targets[i],
                                      crossValidation = crossValidation,
                                      penalty = penalty,
                                      algorithm = algorithm,
                                      nFolds = nFolds,
                                      maxSuppSize = maxSuppSize,
                                      verbose = verbose)
      weightDT <- rbind(weightDT, subWeightDT)
    }

  } else {
    cores <- min(parallel::detectCores(logical = FALSE), cores, length(targets))
    cl <- snow::makeSOCKcluster(cores)
    doSNOW::registerDoSNOW(cl)
    pb <- txtProgressBar(
      min = 1,
      max = length(targets),
      width = 100,
      style = 3
    )
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    "%dopar%" <- foreach::"%dopar%"
    weightDT <- foreach::foreach(
      target = targets,
      .combine = "rbind",
      .export = c("sparse.regression", "sub.inferCSN"),
      .packages = "Kendall",
      .options.snow = opts
    ) %dopar% {
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

    close(pb)
    parallel::stopCluster(cl)
  }

  attr(weightDT, "rng") <- NULL
  attr(weightDT, "doRNG_version") <- NULL
  weightDT <- weightDT[order(abs(weightDT$weight), decreasing = TRUE), ]
  return(weightDT)
}
