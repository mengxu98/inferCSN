#' @title Inferring Cell-Specific Gene Regulatory Network
#'
#' @useDynLib inferCSN
#'
#' @param object Input object
#' @param ... Arguments for other methods
#'
#' @docType methods
#' @rdname inferCSN
#' @return A data table of gene-gene regulatory relationship
#' @export
#'
#' @examples
#' library(inferCSN)
#' data("example_matrix")
#' weight_table <- inferCSN(example_matrix, verbose = TRUE)
#' head(weight_table)
#'
#' weight_table <- inferCSN(example_matrix, verbose = TRUE, cores = 2)
#' head(weight_table)
setGeneric(
  "inferCSN",
  signature = "object",
  function(object, ...) {
    UseMethod(generic = "inferCSN", object = object)
  }
)

#' @param penalty The type of regularization.
#' This can take either one of the following choices: "L0" and "L0L2".
#' For high-dimensional and sparse data, such as single-cell sequencing data, "L0L2" is more effective.
#' @param algorithm The type of algorithm used to minimize the objective function.
#' Currently "CD" and "CDPSI" are supported.
#' The CDPSI algorithm may yield better results, but it also increases running time.
#' @param cross_validation Check whether cross validation is used.
#' @param n_folds The number of folds for cross-validation.
#' @param seed The seed used in randomly shuffling the data for cross-validation.
#' @param k_folds The number of folds for sample split.
#' @param r_threshold Threshold of R^2.
#' @param regulators Regulator genes.
#' @param targets Target genes.
#' @param regulators_num The number of non-zore coef, this value will affect the final performance.
#' The maximum support size at which to terminate the regularization path.
#' Recommend setting this to a small fraction of min(n,p) (e.g. 0.05 * min(n,p)) as L0 regularization typically selects a small portion of non-zeros.
#' @param verbose Print detailed information.
#' @param cores CPU cores.
#'
#' @import Matrix
#' @import doParallel
#' @import parallel
#' @import foreach
#' @import progress
#'
#' @importFrom methods as is
#' @importFrom Rcpp evalCpp
#' @importFrom stats coef predict
#' @importFrom utils methods
#'
#' @rdname inferCSN
#' @export
setMethod(
  "inferCSN",
  signature = "matrix",
  function(
    object,
    penalty = "L0",
    algorithm = "CD",
    cross_validation = FALSE,
    seed = 1,
    n_folds = 10,
    k_folds = NULL,
    r_threshold = 0,
    regulators = NULL,
    targets = NULL,
    regulators_num = NULL,
    verbose = FALSE,
    cores = 1,
    ...) {
    if (verbose) message(paste("Running start for <", class(object)[1], ">."))
    matrix <- object
    rm(object)

    # Check input parameters
    check.parameters(
      matrix = matrix,
      penalty = penalty,
      algorithm = algorithm,
      cross_validation = cross_validation,
      seed = seed,
      n_folds = n_folds,
      k_folds = k_folds,
      r_threshold = r_threshold,
      regulators = regulators,
      targets = targets,
      regulators_num = regulators_num,
      verbose = verbose,
      cores = cores
    )

    if (!is.null(regulators)) {
      regulators_matrix <- matrix[, intersect(colnames(matrix), regulators)]
    } else {
      regulators_matrix <- matrix
    }

    if (!is.null(targets)) {
      targets_matrix <- matrix[, intersect(colnames(matrix), targets)]
    } else {
      targets_matrix <- matrix
    }
    targets <- colnames(targets_matrix)
    rm(matrix)

    target <- NULL
    cores <- min(
      (parallel::detectCores(logical = FALSE) - 1), cores, length(targets)
    )
    if (cores == 1) {
      if (verbose) message("Using 1 core.")
      # Format progress information
      format <- "Running [:bar] :percent, No.:current of :total genes, :elapsed."
      pb <- progress::progress_bar$new(
        format = format,
        total = length(targets),
        clear = TRUE,
        width = 80
      )

      weight_table <- purrr::map_dfr(
        targets, function(target) {
          if (verbose) pb$tick()
          sub.inferCSN(
            regulators_matrix = regulators_matrix,
            targets_matrix = targets_matrix,
            target = target,
            cross_validation = cross_validation,
            seed = seed,
            penalty = penalty,
            algorithm = algorithm,
            n_folds = n_folds,
            k_folds = k_folds,
            r_threshold = r_threshold,
            regulators_num = regulators_num,
            verbose = verbose
          )
        }
      )
    } else {
      doParallel::registerDoParallel(cores = cores)
      if (verbose) message("Using ", foreach::getDoParWorkers(), " cores.")

      "%dopar%" <- foreach::"%dopar%"
      weight_table <- foreach::foreach(
        target = targets,
        .export = c("sub.inferCSN", "sparse.regression")
      ) %dopar% {
        sub.inferCSN(
          regulators_matrix = regulators_matrix,
          targets_matrix = targets_matrix,
          target = target,
          cross_validation = cross_validation,
          seed = seed,
          penalty = penalty,
          algorithm = algorithm,
          n_folds = n_folds,
          k_folds = k_folds,
          r_threshold = r_threshold,
          regulators_num = regulators_num,
          verbose = verbose
        )
      }
      weight_table <- purrr::list_rbind(weight_table)
      doParallel::stopImplicitCluster()
    }

    weight_table <- weight_table[order(
      abs(as.numeric(weight_table$weight)),
      decreasing = TRUE
    ), ]
    if (verbose) message("Run done.")
    return(weight_table)
  }
)

#' @rdname inferCSN
#' @export
setMethod(
  "inferCSN",
  signature = "data.frame",
  function(
    object,
    penalty = "L0",
    algorithm = "CD",
    cross_validation = FALSE,
    seed = 1,
    n_folds = 10,
    k_folds = NULL,
    r_threshold = 0,
    regulators = NULL,
    targets = NULL,
    regulators_num = NULL,
    verbose = FALSE,
    cores = 1,
    ...) {
    if (verbose) warning("Converting class type of input data from <data.frame> to <matrix>.")
    matrix <- as.matrix(object)
    rm(object)

    inferCSN(
      object = matrix,
      penalty = penalty,
      algorithm = algorithm,
      cross_validation = cross_validation,
      seed = seed,
      n_folds = n_folds,
      k_folds = k_folds,
      r_threshold = r_threshold,
      regulators = regulators,
      targets = targets,
      regulators_num = regulators_num,
      verbose = verbose,
      cores = cores
    )
  }
)
