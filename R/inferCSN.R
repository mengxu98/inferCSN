#' @title **infer**ring **C**ell-**S**pecific gene regulatory **N**etwork
#'
#' @useDynLib inferCSN
#'
#' @md
#' @param object The input data for *`inferCSN`*.
#' @param penalty The type of regularization, default is *`L0`*.
#' This can take either one of the following choices: *`L0`*, *`L0L1`*, and *`L0L2`*.
#' For high-dimensional and sparse data, *`L0L2`* is more effective.
#' @param cross_validation Logical value, default is *`FALSE`*, whether to use cross-validation.
#' @param n_folds The number of folds for cross-validation, default is *`5`*.
#' @param seed The random seed for cross-validation, default is *`1`*.
#' @param subsampling_method The method to use for subsampling. Options are "sample", "pseudobulk" or "meta_cells".
#' @param subsampling_ratio The percent of all samples used for \code{\link{fit_srm}}, default is *`1`*.
#' @param r_squared_threshold Threshold of \eqn{R^2} coefficient, default is *`0`*.
#' @param regulators The regulator genes for which to infer the regulatory network.
#' @param targets The target genes for which to infer the regulatory network.
#' Recommend setting this to a small fraction of min(n,p) (e.g. 0.05 * min(n,p)) as L0 regularization typically selects a small portion of non-zeros.
#' @param cores The number of cores to use for parallelization with \code{\link[foreach]{foreach}}, default is *`1`*.
#' @param verbose Logical value, default is *`TRUE`*, whether to print progress messages.
#' @param ... Parameters for other methods.
#'
#' @docType methods
#' @rdname inferCSN
#' @return A data table of regulator-target regulatory relationships
#' @export
setGeneric(
  name = "inferCSN",
  signature = c("object"),
  def = function(object,
                 penalty = "L0",
                 cross_validation = FALSE,
                 seed = 1,
                 n_folds = 5,
                 subsampling_method = c(
                   "sample", "meta_cells", "pseudobulk"
                 ),
                 subsampling_ratio = 1,
                 r_squared_threshold = 0,
                 regulators = NULL,
                 targets = NULL,
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
#'
#' network_table_2 <- inferCSN(
#'   example_matrix,
#'   cores = 2
#' )
#'
#' head(network_table_1)
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
#'
#' \dontrun{
#' data("example_ground_truth")
#' network_table_07 <- inferCSN(
#'   example_matrix,
#'   r_squared_threshold = 0.7
#' )
#' calculate_metrics(
#'   network_table_1,
#'   example_ground_truth,
#'   plot = TRUE
#' )
#' calculate_metrics(
#'   network_table_07,
#'   example_ground_truth,
#'   plot = TRUE
#' )
#' }
setMethod(
  f = "inferCSN",
  signature = signature(object = "matrix"),
  definition = function(object,
                        penalty = "L0",
                        cross_validation = FALSE,
                        seed = 1,
                        n_folds = 5,
                        subsampling_method = c(
                          "sample", "meta_cells", "pseudobulk"
                        ),
                        subsampling_ratio = 1,
                        r_squared_threshold = 0,
                        regulators = NULL,
                        targets = NULL,
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
      cross_validation = cross_validation,
      seed = seed,
      n_folds = n_folds,
      subsampling_method = subsampling_method,
      subsampling_ratio = subsampling_ratio,
      r_squared_threshold = r_squared_threshold,
      regulators = regulators,
      targets = targets,
      verbose = verbose,
      cores = cores,
      ...
    )

    object <- subsampling(
      matrix = object,
      subsampling_method = subsampling_method,
      subsampling_ratio = subsampling_ratio,
      seed = seed,
      verbose = verbose,
      ...
    )

    regulators <- intersect(
      colnames(object),
      regulators %ss% colnames(object)
    )
    targets <- intersect(
      colnames(object),
      targets %ss% colnames(object)
    )

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
          n_folds = n_folds,
          subsampling_ratio = subsampling_ratio,
          r_squared_threshold = r_squared_threshold,
          verbose = verbose,
          ...
        )
      },
      cores = cores,
      verbose = verbose
    ) |>
      purrr::list_rbind() |>
      network_format(abs_weight = FALSE)

    log_message(
      "Run done.",
      message_type = "success",
      verbose = verbose
    )

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
                        cross_validation = FALSE,
                        seed = 1,
                        n_folds = 5,
                        subsampling_method = c(
                          "sample", "meta_cells", "pseudobulk"
                        ),
                        subsampling_ratio = 1,
                        r_squared_threshold = 0,
                        regulators = NULL,
                        targets = NULL,
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
      cross_validation = cross_validation,
      seed = seed,
      n_folds = n_folds,
      subsampling_method = subsampling_method,
      subsampling_ratio = subsampling_ratio,
      r_squared_threshold = r_squared_threshold,
      regulators = regulators,
      targets = targets,
      verbose = verbose,
      cores = cores,
      ...
    )

    object <- subsampling(
      matrix = object,
      subsampling_method = subsampling_method,
      subsampling_ratio = subsampling_ratio,
      seed = seed,
      verbose = verbose,
      ...
    )

    regulators <- intersect(
      colnames(object),
      regulators %ss% colnames(object)
    )
    targets <- intersect(
      colnames(object),
      targets %ss% colnames(object)
    )

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
          n_folds = n_folds,
          r_squared_threshold = r_squared_threshold,
          verbose = verbose,
          ...
        )
      },
      cores = cores,
      verbose = verbose
    ) |>
      purrr::list_rbind() |>
      network_format(abs_weight = FALSE)

    log_message(
      "Run done.",
      message_type = "success",
      verbose = verbose
    )

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
                        cross_validation = FALSE,
                        seed = 1,
                        n_folds = 5,
                        subsampling_method = c(
                          "sample", "meta_cells", "pseudobulk"
                        ),
                        subsampling_ratio = 1,
                        r_squared_threshold = 0,
                        regulators = NULL,
                        targets = NULL,
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
      cross_validation = cross_validation,
      seed = seed,
      n_folds = n_folds,
      subsampling_method = subsampling_method,
      subsampling_ratio = subsampling_ratio,
      r_squared_threshold = r_squared_threshold,
      regulators = regulators,
      targets = targets,
      verbose = verbose,
      cores = cores,
      ...
    )
  }
)
