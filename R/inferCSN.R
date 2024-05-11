#' @title Inferring Cell-Specific Gene Regulatory Network
#'
#' @useDynLib inferCSN
#'
#' @param object Input object
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
#' @param cores CPU cores.
#' @param verbose Print detailed information.
#' @param ... Arguments for other methods
#'
#' @docType methods
#' @rdname inferCSN
#' @return A data table of gene-gene regulatory relationship
#' @export
setGeneric(
  name = "inferCSN",
  signature = c("object"),
  def = function(
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
    cores = 1,
    verbose = FALSE,
    ...) {
    UseMethod(
      generic = "inferCSN",
      object = object
    )
  }
)

#' @rdname inferCSN
#' @export
#'
#' @examples
#' data("example_matrix")
#' weight_table <- inferCSN(example_matrix, verbose = TRUE)
#' head(weight_table)
#'
#' weight_table <- inferCSN(example_matrix, cores = 2)
#' head(weight_table)
setMethod(
  f = "inferCSN",
  signature = signature(
    object = "matrix"),
  definition = function(
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
    cores = 1,
    verbose = FALSE,
    ...) {
    if (verbose) {
      message(paste("Running start for <", class(object)[1], ">."))
    }

    # Check input parameters
    check.parameters(
      matrix = object,
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
      regulators <- intersect(colnames(object), regulators)
    } else {
      regulators <- colnames(object)
    }

    if (!is.null(targets)) {
      targets <- intersect(colnames(object), targets)
    } else {
      targets <- colnames(object)
    }

    names(targets) <- targets
    cores <- min(
      (parallel::detectCores(logical = FALSE) - 1), cores, length(targets)
    )
    weight_list <- parallelize_fun(
      x = targets,
      fun = function(x) {
        single.network(
          matrix = object,
          regulators = regulators,
          target = x,
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
      },
      cores = cores,
      verbose = verbose
    )
    weight_table <- purrr::list_rbind(weight_list)
    weight_table <- net.format(
      weight_table,
      abs_weight = FALSE
    )
    if (verbose) message("Run done.")

    return(weight_table)
  }
)

#' @rdname inferCSN
#' @export
#'
#' @examples
#' data("example_matrix")
#' weight_table <- inferCSN(
#'   as.data.frame(example_matrix),
#'   verbose = TRUE
#' )
#' head(weight_table)
setMethod(
  f = "inferCSN",
  signature = signature(
    object = "data.frame"
  ),
  definition = function(
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
    cores = 1,
    verbose = FALSE,
    ...) {
    if (verbose) {
      warning("Converting class type of input data from <data.frame> to <matrix>.")
    }
    object <- as.matrix(object)

    inferCSN(
      object = object,
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
      cores = cores,
      ...
    )
  }
)
