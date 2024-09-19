#' @title **infer**ring **C**ell-**S**pecific gene regulatory **N**etwork
#'
#' @useDynLib inferCSN
#'
#' @param object The input data for *`inferCSN`*.
#' @param penalty The type of regularization, default is *`L0`*.
#' This can take either one of the following choices: *`L0`*, *`L0L1`*, and *`L0L2`*.
#' For high-dimensional and sparse data, *`L0L2`* is more effective.
#' @param algorithm The type of algorithm used to minimize the objective function, default is *`CD`*.
#' Currently *`CD`* and *`CDPSI`* are supported.
#' The *`CDPSI`* algorithm may yield better results, but it also increases running time.
#' @param cross_validation Logical value, default is *`FALSE`*, whether to use cross-validation.
#' @param n_folds The number of folds for cross-validation, default is *`10`*.
#' @param seed The random seed for cross-validation, default is *`1`*.
#' @param subsampling The percent of all samples used for \code{\link{sparse_regression}}, default is *`1`*.
#' @param r_threshold Threshold of \eqn{R^2} or correlation coefficient, default is *`0`*.
#' @param regulators The regulator genes for which to infer the regulatory network.
#' @param targets The target genes for which to infer the regulatory network.
#' @param regulators_num The number of non-zore coefficients, this value will affect the final performance.
#' The maximum support size at which to terminate the regularization path.
#' Recommend setting this to a small fraction of min(n,p) (e.g. 0.05 * min(n,p)) as L0 regularization typically selects a small portion of non-zeros.
#' @param cores The number of cores to use for parallelization with \code{\link[foreach]{foreach}}, default is *`1`*.
#' @param verbose Logical value, default is *`TRUE`*, whether to print progress messages.
#' @param ... Parameters for other methods.
#'
#' @md
#' @docType methods
#' @rdname inferCSN
#' @return A data table of regulator-target regulatory relationships
#' @export
setGeneric(
  name = "inferCSN",
  signature = c("object"),
  def = function(object,
                 penalty = "L0",
                 algorithm = "CD",
                 cross_validation = FALSE,
                 seed = 1,
                 n_folds = 10,
                 subsampling = 1,
                 r_threshold = 0,
                 regulators = NULL,
                 targets = NULL,
                 regulators_num = NULL,
                 cores = 1,
                 verbose = TRUE,
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
#' network_table_1 <- inferCSN(
#'   example_matrix
#' )
#' head(network_table_1)
#'
#' network_table_2 <- inferCSN(
#'   example_matrix,
#'   cores = 2
#' )
#'
#' identical(
#'   network_table_1,
#'   network_table_2
#' )
#'
#' inferCSN(
#'   example_matrix,
#'   regulators = c("g1", "g2"),
#'   targets = c("g3", "g4")
#' )
#' inferCSN(
#'   example_matrix,
#'   regulators = c("g1", "g2"),
#'   targets = c("g3", "g0")
#' )
setMethod(
  f = "inferCSN",
  signature = signature(object = "matrix"),
  definition = function(object,
                        penalty = "L0",
                        algorithm = "CD",
                        cross_validation = FALSE,
                        seed = 1,
                        n_folds = 10,
                        subsampling = 1,
                        r_threshold = 0,
                        regulators = NULL,
                        targets = NULL,
                        regulators_num = NULL,
                        cores = 1,
                        verbose = TRUE,
                        ...) {
    log_message(
      "Running for <dense matrix>.",
      verbose = verbose
    )

    .check_parameters(
      matrix = object,
      penalty = penalty,
      algorithm = algorithm,
      cross_validation = cross_validation,
      seed = seed,
      n_folds = n_folds,
      subsampling = subsampling,
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

    network_table <- parallelize_fun(
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
          subsampling = subsampling,
          r_threshold = r_threshold,
          regulators_num = regulators_num,
          verbose = verbose,
          ...
        )
      },
      cores = cores,
      verbose = verbose
    ) |>
      purrr::list_rbind() |>
      network_format(abs_weight = FALSE)

    log_message("Run done.", verbose = verbose)

    return(network_table)
  }
)

#' @rdname inferCSN
#' @export
#'
#' @examples
#' \dontrun{
#' data("example_matrix")
#' network_table <- inferCSN(example_matrix)
#' head(network_table)
#'
#' network_table_sparse_1 <- inferCSN(
#'   as(example_matrix, "sparseMatrix")
#' )
#' head(network_table_sparse_1)
#'
#' network_table_sparse_2 <- inferCSN(
#'   as(example_matrix, "sparseMatrix"),
#'   cores = 2
#' )
#' identical(
#'   network_table,
#'   network_table_sparse_1
#' )
#'
#' identical(
#'   network_table_sparse_1,
#'   network_table_sparse_2
#' )
#'
#' plot_scatter(
#'   data.frame(
#'     network_table$weight,
#'     network_table_sparse_1$weight
#'   ),
#'   legend_position = "none"
#' )
#'
#' plot_weight_distribution(
#'   network_table
#' ) + plot_weight_distribution(
#'   network_table_sparse_1
#' )
#' }
setMethod(
  f = "inferCSN",
  signature = signature(object = "sparseMatrix"),
  definition = function(object,
                        penalty = "L0",
                        algorithm = "CD",
                        cross_validation = FALSE,
                        seed = 1,
                        n_folds = 10,
                        subsampling = 1,
                        r_threshold = 0,
                        regulators = NULL,
                        targets = NULL,
                        regulators_num = NULL,
                        cores = 1,
                        verbose = TRUE,
                        ...) {
    log_message(
      "Running for <", class(object), ">.",
      verbose = verbose
    )

    .check_parameters(
      matrix = object,
      penalty = penalty,
      algorithm = algorithm,
      cross_validation = cross_validation,
      seed = seed,
      n_folds = n_folds,
      subsampling = subsampling,
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

    network_table <- parallelize_fun(
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
          subsampling = subsampling,
          r_threshold = r_threshold,
          regulators_num = regulators_num,
          verbose = verbose,
          ...
        )
      },
      cores = cores,
      verbose = verbose
    ) |>
      purrr::list_rbind() |>
      network_format(abs_weight = FALSE)

    log_message("Run done.", verbose = verbose)

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
                        subsampling = 1,
                        r_threshold = 0,
                        regulators = NULL,
                        targets = NULL,
                        regulators_num = NULL,
                        cores = 1,
                        verbose = TRUE,
                        ...) {
    log_message(
      "convert the class type of the input data from <data.frame> to <matrix>.",
      message_type = "warning",
      verbose = verbose
    )

    inferCSN(
      object = as_matrix(object),
      penalty = penalty,
      algorithm = algorithm,
      cross_validation = cross_validation,
      seed = seed,
      n_folds = n_folds,
      subsampling = subsampling,
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
