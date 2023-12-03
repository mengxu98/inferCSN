#' @useDynLib inferCSN

#' @title Inferring Cell-Specific Gene Regulatory Network
#'
#' @param matrix An expression matrix, cells by genes
#' @param penalty The type of regularization.
#' This can take either one of the following choices: "L0" and "L0L2".
#' For high-dimensional and sparse data, such as single-cell sequencing data, "L0L2" is more effective.
#' @param algorithm The type of algorithm used to minimize the objective function.
#' Currently "CD" and "CDPSI" are supported.
#' The CDPSI algorithm may yield better results, but it also increases running time.
#' @param crossValidation Check whether cross validation is used.
#' @param nFolds The number of folds for cross-validation.
#' @param seed The seed used in randomly shuffling the data for cross-validation.
#' @param kFolds The number of folds for sample split.
#' @param rThreshold Threshold of R^2.
#' @param regulators Regulator genes.
#' @param targets Target genes.
#' @param maxSuppSize The number of non-zore coef, this value will affect the final performance.
#' The maximum support size at which to terminate the regularization path.
#' Recommend setting this to a small fraction of min(n,p) (e.g. 0.05 * min(n,p)) as L0 regularization typically selects a small portion of non-zeros.
#' @param verbose Print detailed information.
#' @param cores CPU cores.
#' @param ... Arguments passed to other methods.
#'
#' @import Matrix
#'
#' @importFrom methods as is
#' @importFrom Rcpp evalCpp
#' @importFrom stats coef predict
#' @importFrom utils methods
#'
#' @return A data table of gene-gene regulatory relationship
#' @export
#'
#' @examples
#' library(inferCSN)
#' data("exampleMatrix")
#' weightDT <- inferCSN(exampleMatrix, verbose = TRUE)
#' head(weightDT)
#'
#' weightDT <- inferCSN(exampleMatrix, verbose = TRUE, cores = 2)
#' head(weightDT)
#'

#' @docType methods
#' @rdname inferCSN
#' @export
#'
setGeneric("inferCSN",
           signature = "matrix",
           function(matrix, ...) {
             UseMethod(generic = "inferCSN", object = matrix)
           })

#' @rdname inferCSN
#' @export
#'
setMethod("inferCSN",
          signature = "data.frame",
          function(matrix,
                   penalty = "L0",
                   algorithm = "CD",
                   crossValidation = FALSE,
                   seed = 1,
                   nFolds = 10,
                   kFolds = NULL,
                   rThreshold = 0,
                   regulators = NULL,
                   targets = NULL,
                   maxSuppSize = NULL,
                   verbose = FALSE,
                   cores = 1,
                   ...) {
            warning("Converting the class type of input data from <data.frame> to <matrix>.")
            matrix <- as.matrix(matrix)
            .inferCSN(matrix = matrix,
                      penalty = penalty,
                      algorithm = algorithm,
                      crossValidation = crossValidation,
                      seed = 1,
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
#' @export
#'
setMethod("inferCSN",
          signature = "matrix",
          function(matrix,
                   penalty = "L0",
                   algorithm = "CD",
                   crossValidation = FALSE,
                   seed = 1,
                   nFolds = 10,
                   kFolds = NULL,
                   rThreshold = 0,
                   regulators = NULL,
                   targets = NULL,
                   maxSuppSize = NULL,
                   verbose = FALSE,
                   cores = 1,
                   ...) {
            .inferCSN(matrix = matrix,
                      penalty = penalty,
                      algorithm = algorithm,
                      crossValidation = crossValidation,
                      seed = seed,
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
                      seed,
                      nFolds,
                      kFolds,
                      rThreshold,
                      regulators,
                      targets,
                      maxSuppSize,
                      verbose,
                      cores,
                      ...) {
  if(verbose) message("Running start.")

  # Check input parameters
  check.parameters(matrix = matrix,
                   penalty = penalty,
                   algorithm = algorithm,
                   crossValidation = crossValidation,
                   seed = seed,
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

  if (!is.null(targets)) {
    targetsMatrix <- matrix[, intersect(colnames(matrix), targets)]
  } else {
    targetsMatrix <- matrix
  }
  targets <- colnames(targetsMatrix)
  rm(matrix)

  cores <- min((parallel::detectCores(logical = FALSE) - 1), cores, length(targets))
  if (cores == 1) {
    if(verbose) message("Using 1 core.")
    # Format progress information
    format <- "Running [:bar] :percent, No.:current of :total genes, :elapsed."
    pb <- progress::progress_bar$new(format = format,
                                     total = length(targets),
                                     clear = TRUE,
                                     width = 80)

    weightDT <- purrr::map_dfr(targets, function(x) {
      if (verbose) pb$tick()
      sub.inferCSN(regulatorsMatrix = regulatorsMatrix,
                   targetsMatrix = targetsMatrix,
                   target = x,
                   crossValidation = crossValidation,
                   seed = seed,
                   penalty = penalty,
                   algorithm = algorithm,
                   nFolds = nFolds,
                   kFolds = kFolds,
                   rThreshold = rThreshold,
                   maxSuppSize = maxSuppSize,
                   verbose = verbose)
    })

  } else {
    doParallel::registerDoParallel(cores = cores)
    if(verbose) message("Using ", foreach::getDoParWorkers(), " cores.")

    "%dopar%" <- foreach::"%dopar%"
    weightDT <- foreach::foreach(target = targets,
                                 .export = c("sub.inferCSN", "sparse.regression")) %dopar% {
                                   sub.inferCSN(regulatorsMatrix = regulatorsMatrix,
                                                targetsMatrix = targetsMatrix,
                                                target = target,
                                                crossValidation = crossValidation,
                                                seed = seed,
                                                penalty = penalty,
                                                algorithm = algorithm,
                                                nFolds = nFolds,
                                                kFolds = kFolds,
                                                rThreshold = rThreshold,
                                                maxSuppSize = maxSuppSize,
                                                verbose = verbose)
                                 }
    weightDT <- data.table::rbindlist(weightDT)
    attr(weightDT, ".internal.selfref") <- NULL
    doParallel::stopImplicitCluster()
  }

  weightDT <- weightDT[order(abs(as.numeric(weightDT$weight)), decreasing = TRUE), ]
  if (verbose) message("Run done.")
  return(weightDT)
}
