#' @title Inferring cell-specific gene regulatory network
#'
#' @param matrix An expression matrix, cells by genes
#' @param penalty The type of regularization.
#' This can take either one of the following choices: "L0"and "L0L2".
#' For high-dimensional and sparse data, such as single-cell transcriptome data, "L0L2" is more effective
#' @param algorithm Currently "CD" and "CDPSI" are supported.
#' The CDPSI algorithm may yield better results, but it also increases running time
#' @param crossValidation Check whether cross validation is used
#' @param nFolds The number of folds for cross-validation
#' @param kFolds The number of folds for sample split
#' @param rThreshold rThreshold
#' @param regulators Regulator genes
#' @param targets Target genes
#' @param maxSuppSize The number of non-zore coef, this value will affect the final performance
#' @param verbose Print detailed information
#' @param cores CPU cores
#'
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

# @docType methods
#' @rdname inferCSN
#' @export
#'
setGeneric("inferCSN",
           signature = "matrix",
           function(matrix,
                    penalty = "L0",
                    algorithm = "CD",
                    crossValidation = FALSE,
                    nFolds = 10,
                    kFolds = NULL,
                    rThreshold = 0,
                    regulators = NULL,
                    targets = NULL,
                    maxSuppSize = NULL,
                    verbose = FALSE,
                    cores = 1) {
             standardGeneric("inferCSN")
           })

#' @rdname inferCSN
#' @aliases inferCSN, data.frame-method
#' @export
#'
setMethod("inferCSN",
          "data.frame",
          function(matrix,
                   penalty = "L0",
                   algorithm = "CD",
                   crossValidation = FALSE,
                   nFolds = 10,
                   kFolds = NULL,
                   rThreshold = 0,
                   regulators = NULL,
                   targets = NULL,
                   maxSuppSize = NULL,
                   verbose = FALSE,
                   cores = 1) {
            matrix <- as.matrix(matrix)
            .inferCSN(matrix = matrix,
                      penalty = penalty,
                      algorithm = algorithm,
                      crossValidation = crossValidation,
                      nFolds = nFolds,
                      kFolds = kFolds,
                      rThreshold = rThreshold,
                      regulators = regulators,
                      targets = targets,
                      maxSuppSize = maxSuppSize,
                      verbose = verbose,
                      cores = cores)
          })

#' @rdname inferCSN
#' @aliases inferCSN, matrix-method
#' @export
#'
setMethod("inferCSN",
          "matrix",
          function(matrix,
                   penalty = "L0",
                   algorithm = "CD",
                   crossValidation = FALSE,
                   nFolds = 10,
                   kFolds = NULL,
                   rThreshold = 0,
                   regulators = NULL,
                   targets = NULL,
                   maxSuppSize = NULL,
                   verbose = FALSE,
                   cores = 1) {
            .inferCSN(matrix = matrix,
                      penalty = penalty,
                      algorithm = algorithm,
                      crossValidation = crossValidation,
                      nFolds = nFolds,
                      kFolds = kFolds,
                      rThreshold = rThreshold,
                      regulators = regulators,
                      targets = targets,
                      maxSuppSize = maxSuppSize,
                      verbose = verbose,
                      cores = cores)
          })

.inferCSN <- function(matrix,
                      penalty,
                      algorithm,
                      crossValidation,
                      nFolds,
                      kFolds,
                      rThreshold,
                      regulators,
                      targets,
                      maxSuppSize,
                      verbose,
                      cores) {
  # Check input arguments
  check.arguments(matrix = matrix,
                  penalty = penalty,
                  algorithm = algorithm,
                  crossValidation = crossValidation,
                  nFolds = nFolds,
                  kFolds = kFolds,
                  rThreshold = rThreshold,
                  regulators = regulators,
                  targets = targets,
                  maxSuppSize = maxSuppSize,
                  verbose = verbose,
                  cores = cores)

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
  rm(matrix)

  if (verbose) message("Using '", penalty, "' penalty......")
  if (verbose & crossValidation) message("Using '", penalty, "' penalty and cross validation......")

  cores <- min(parallel::detectCores(logical = FALSE), cores, length(targets))
  if (cores == 1) {
    # Format progress information
    pb <- progress::progress_bar$new(format = "Running [:bar] :percent, No.:current of :total gene,:elapsed......",
                                     total = length(targets),
                                     clear = TRUE,
                                     width = 100)

    weightDT <- purrr::map_dfr(regulators, function(x) {
      if (verbose) pb$tick()
      sub.inferCSN(regulatorsMatrix = regulatorsMatrix,
                   targetsMatrix = targetsMatrix,
                   target = x,
                   crossValidation = crossValidation,
                   penalty = penalty,
                   algorithm = algorithm,
                   nFolds = nFolds,
                   kFolds = kFolds,
                   rThreshold = rThreshold,
                   maxSuppSize = maxSuppSize,
                   verbose = verbose)
    })

  } else {
    cl <- parallel::makeCluster(cores)
    doSNOW::registerDoSNOW(cl)

    if (verbose) {
      pb <- txtProgressBar(min = 1,
                           max = length(targets),
                           width = 100,
                           style = 3)

      progress <- function(n) setTxtProgressBar(pb, n)
      opts <- list(progress = progress)
    } else {
      opts <- NULL
    }

    "%dopar%" <- foreach::"%dopar%"
    weightDT <- foreach::foreach(target = targets,
                                 .export = c("sub.inferCSN", "sparse.regression"),
                                 .options.snow = opts) %dopar% {
                                   sub.inferCSN(regulatorsMatrix = regulatorsMatrix,
                                                targetsMatrix = targetsMatrix,
                                                target = target,
                                                crossValidation = crossValidation,
                                                penalty = penalty,
                                                algorithm = algorithm,
                                                nFolds = nFolds,
                                                kFolds = kFolds,
                                                rThreshold = rThreshold,
                                                maxSuppSize = maxSuppSize,
                                                verbose = verbose)
                                 }
    weightDT <- purrr::list_rbind(weightDT)

    if (verbose) close(pb)
    doParallel::stopImplicitCluster()
    parallel::stopCluster(cl)
  }

  weightDT <- weightDT[order(abs(as.numeric(weightDT$weight)), decreasing = TRUE), ]
  return(weightDT)
}
