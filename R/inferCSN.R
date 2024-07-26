#' @title Inferring Cell-Specific Gene Regulatory Network
#'
#' @useDynLib inferCSN
#'
#' @param object The input data for \code{inferCSN}.
#' @param penalty The type of regularization.
#' This can take either one of the following choices: \code{L0}, \code{L0L1} and \code{L0L2}.
#' For high-dimensional and sparse data, such as single-cell sequencing data, \code{L0L2} is more effective.
#' @param algorithm The type of algorithm used to minimize the objective function.
#' Currently \code{CD} and \code{CDPSI} are supported.
#' The \code{CDPSI} algorithm may yield better results, but it also increases running time.
#' @param cross_validation Check whether cross validation is used.
#' @param n_folds The number of folds for cross-validation.
#' @param seed The seed used in randomly shuffling the data for cross-validation.
#' @param percent_samples The percent of all samples used for \code{\link{sparse_regression}}. Default set to 1.
#' @param r_threshold Threshold of \eqn{R^2} or correlation coefficient.
#' @param regulators A character vector with the regulators to consider for CSN inference.
#' @param targets A character vector with the targets to consider for CSN inference.
#' @param regulators_num The number of non-zore coefficients, this value will affect the final performance.
#' The maximum support size at which to terminate the regularization path.
#' Recommend setting this to a small fraction of min(n,p) (e.g. 0.05 * min(n,p)) as L0 regularization typically selects a small portion of non-zeros.
#' @param cores Number of CPU cores used. Setting to parallelize the computation with \code{\link[foreach]{foreach}}.
#' @param verbose Logical value. Whether to print detailed information.
#' @param ... Parameters for other methods.
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
      percent_samples = 1,
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
#' network_table <- inferCSN(example_matrix, verbose = TRUE)
#' head(network_table)
#'
#' network_table <- inferCSN(example_matrix, cores = 2)
#' head(network_table)
setMethod(
  f = "inferCSN",
  signature = signature(object = "matrix"),
  definition = function(object,
                        penalty = "L0",
                        algorithm = "CD",
                        cross_validation = FALSE,
                        seed = 1,
                        n_folds = 10,
                        percent_samples = 1,
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
    check_parameters(
      matrix = object,
      penalty = penalty,
      algorithm = algorithm,
      cross_validation = cross_validation,
      seed = seed,
      n_folds = n_folds,
      percent_samples = percent_samples,
      r_threshold = r_threshold,
      regulators = regulators,
      targets = targets,
      regulators_num = regulators_num,
      verbose = verbose,
      cores = cores,
      ...
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
    if (is.null(regulators_num)) {
      regulators_num <- (ncol(object) - 1)
    }
    names(targets) <- targets
    cores <- .cores_detect(cores, length(targets))

    weight_list <- parallelize_fun(
      x = targets,
      fun = function(x) {
        single_network(
          matrix = object,
          regulators = regulators,
          target = x,
          cross_validation = cross_validation,
          seed = seed,
          penalty = penalty,
          algorithm = algorithm,
          n_folds = n_folds,
          percent_samples = percent_samples,
          r_threshold = r_threshold,
          regulators_num = regulators_num,
          verbose = verbose
        )
      },
      cores = cores,
      verbose = verbose
    )
    network_table <- purrr::list_rbind(weight_list)
    network_table <- network_format(
      network_table,
      abs_weight = FALSE
    )
    if (verbose) message("Run done.")

    return(network_table)
  }
)

#' @rdname inferCSN
#' @export
#'
#' @examples
#' data("example_matrix")
#' network_table <- inferCSN(example_matrix, verbose = TRUE)
#' head(network_table)
#'
#' network_table <- inferCSN(example_matrix, cores = 2)
#' head(network_table)
#'
#' network_table_sparse <- inferCSN(
#'   as(example_matrix, "sparseMatrix"),
#'   cores = 2
#' )
#' head(network_table_sparse)
#'
#' plot_scatter(
#'   data.frame(
#'     network_table_sparse$weight,
#'     network_table$weight
#'   )
#' )
setMethod(
  f = "inferCSN",
  signature = signature(object = "sparseMatrix"),
  definition = function(object,
                        penalty = "L0",
                        algorithm = "CD",
                        cross_validation = FALSE,
                        seed = 1,
                        n_folds = 10,
                        percent_samples = 1,
                        r_threshold = 0,
                        regulators = NULL,
                        targets = NULL,
                        regulators_num = NULL,
                        cores = 1,
                        verbose = FALSE,
                        ...) {
    if (verbose) {
      message(
        paste0(
          "The class type of input data is <", class(object), ">."
        )
      )
    }

    check_parameters(
      matrix = object,
      penalty = penalty,
      algorithm = algorithm,
      cross_validation = cross_validation,
      seed = seed,
      n_folds = n_folds,
      percent_samples = percent_samples,
      r_threshold = r_threshold,
      regulators = regulators,
      targets = targets,
      regulators_num = regulators_num,
      verbose = verbose,
      cores = cores,
      ...
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
    if (is.null(regulators_num)) {
      regulators_num <- (ncol(object) - 1)
    }
    names(targets) <- targets
    cores <- .cores_detect(cores, length(targets))

    weight_list <- parallelize_fun(
      x = targets,
      fun = function(x) {
        single_network(
          matrix = object,
          regulators = regulators,
          target = x,
          cross_validation = cross_validation,
          seed = seed,
          penalty = penalty,
          algorithm = algorithm,
          n_folds = n_folds,
          percent_samples = percent_samples,
          r_threshold = r_threshold,
          regulators_num = regulators_num,
          verbose = verbose
        )
      },
      cores = cores,
      verbose = verbose
    )
    network_table <- purrr::list_rbind(weight_list)
    network_table <- network_format(
      network_table,
      abs_weight = FALSE
    )
    if (verbose) message("Run done.")

    return(network_table)
  }
)

#' @rdname inferCSN
#' @export
setMethod(
  f = "inferCSN",
  signature = signature(object = "data.frame"),
  definition = function(object,
                        penalty = "L0",
                        algorithm = "CD",
                        cross_validation = FALSE,
                        seed = 1,
                        n_folds = 10,
                        percent_samples = 1,
                        r_threshold = 0,
                        regulators = NULL,
                        targets = NULL,
                        regulators_num = NULL,
                        cores = 1,
                        verbose = FALSE,
                        ...) {
    if (verbose) {
      message(
        paste0(
          "Converting class type of input data from <data.frame> to <matrix>."
        )
      )
    }

    inferCSN(
      object = as_matrix(object),
      penalty = penalty,
      algorithm = algorithm,
      cross_validation = cross_validation,
      seed = seed,
      n_folds = n_folds,
      percent_samples = percent_samples,
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
