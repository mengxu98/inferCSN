#' @title inferring cell-type specific gene regulatory network
#'
#' @md
#' @param object The input data for inferring network.
#' @param penalty The type of regularization.
#' This can take either one of the following choices: `"L0"`, `"L0L1"`, and `"L0L2"`.
#' For high-dimensional and sparse data, `"L0L2"` is more effective.
#' Default is `"L0"`.
#' @param cross_validation Whether to use cross-validation.
#' Default is `FALSE`.
#' @param n_folds The number of folds for cross-validation.
#' Default is `5`.
#' @param seed The random seed for cross-validation.
#' Default is `1`.
#' @param subsampling_method The method to use for subsampling.
#' Options are `"sample"`, `"pseudobulk"` or `"meta_cells"`.
#' @param subsampling_ratio The percent of all samples used for [fit_srm].
#' Default is `1`.
#' @param r_squared_threshold Threshold of R² coefficient.
#' Default is `0`.
#' @param regulators The regulator genes for which to infer the regulatory network.
#' @param targets The target genes for which to infer the regulatory network.
#' @param cores The number of cores to use for parallelization with [foreach::foreach].
#' Default is `1`.
#' @param verbose Whether to print progress messages.
#' Default is `TRUE`.
#' @param ... Parameters for other methods.
#'
#' @docType methods
#' @rdname inferCSN
#' @return
#' A data table of regulator-target regulatory relationships,
#' which has three columns: `regulator`, `target`, and `weight`.
#'
#' @export
setGeneric(
  name = "inferCSN",
  signature = c("object"),
  def = function(object,
                 penalty = c("L0", "L0L1", "L0L2"),
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
#' data(example_matrix)
#' network_table <- inferCSN(example_matrix)
#'
#' head(network_table)
#'
#' inferCSN(
#'   example_matrix,
#'   regulators = c("g1", "g2"),
#'   targets = c("g3", "g4")
#' )
#'
#' inferCSN(
#'   example_matrix,
#'   regulators = c("g1", "g2"),
#'   targets = c("g3", "g0")
#' )
setMethod(
  f = "inferCSN",
  signature = signature(object = "matrix"),
  definition = function(object,
                        penalty = c("L0", "L0L1", "L0L2"),
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
    thisutils::log_message(
      "Inferring network for {.cls {class(object)}}...",
      verbose = verbose
    )

    validated_params <- do.call(
      check_parameters,
      c(
        list(matrix = object),
        mget(names(formals())[-1], environment())
      )
    )

    object <- subsampling(
      matrix = object,
      subsampling_method = validated_params$subsampling_method,
      subsampling_ratio = validated_params$subsampling_ratio,
      seed = validated_params$seed,
      verbose = verbose,
      ...
    )

    if (!is.null(validated_params$regulators)) {
      regulators <- validated_params$regulators
    } else {
      regulators <- regulators %ss% colnames(object)
    }

    if (!is.null(validated_params$targets)) {
      targets <- validated_params$targets
    } else {
      targets <- targets %ss% colnames(object)
    }

    names(targets) <- targets
    param_names <- c(
      "cross_validation", "seed", "penalty",
      "n_folds", "r_squared_threshold", "verbose"
    )
    network_params <- validated_params[param_names]
    param_index <- !vapply(network_params, is.null, logical(1))
    network_params <- network_params[param_index]

    dots_list <- list(...)
    dots_list <- dots_list[!vapply(dots_list, is.null, logical(1))]

    single_network_args <- c(network_params, dots_list)

    network_table <- thisutils::parallelize_fun(
      x = targets,
      fun = function(x) {
        thisutils::invoke_fun(
          single_network,
          c(
            list(
              matrix = object,
              regulators = regulators,
              target = x
            ),
            single_network_args
          )
        )
      },
      cores = validated_params$cores,
      verbose = verbose
    ) |>
      purrr::list_rbind() |>
      network_format(abs_weight = FALSE)

    thisutils::log_message(
      "Inferring network done",
      message_type = "success",
      verbose = verbose
    )
    network_info <- data.frame(
      Edges = nrow(network_table),
      Regulators = length(unique(network_table$regulator)),
      Targets = length(unique(network_table$target))
    )
    thisutils::log_message(
      "Network information:\n",
      network_info,
      text_color = "grey",
      timestamp_style = FALSE,
      verbose = verbose
    )

    return(network_table)
  }
)

#' @rdname inferCSN
#' @export
setMethod(
  f = "inferCSN",
  signature = signature(object = "sparseMatrix"),
  definition = function(object,
                        penalty = c("L0", "L0L1", "L0L2"),
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
    thisutils::log_message(
      "Inferring network for {.cls {class(object)}}...",
      verbose = verbose
    )

    validated_params <- do.call(
      check_parameters,
      c(
        list(matrix = object),
        mget(names(formals())[-1], environment())
      )
    )

    object <- subsampling(
      matrix = object,
      subsampling_method = validated_params$subsampling_method,
      subsampling_ratio = validated_params$subsampling_ratio,
      seed = validated_params$seed,
      verbose = verbose,
      ...
    )

    if (!is.null(validated_params$regulators)) {
      regulators <- validated_params$regulators
    } else {
      regulators <- regulators %ss% colnames(object)
    }

    if (!is.null(validated_params$targets)) {
      targets <- validated_params$targets
    } else {
      targets <- targets %ss% colnames(object)
    }

    names(targets) <- targets

    param_names <- c(
      "cross_validation", "seed", "penalty",
      "n_folds", "r_squared_threshold", "verbose"
    )
    network_params <- validated_params[param_names]
    param_index <- !vapply(network_params, is.null, logical(1))
    network_params <- network_params[param_index]

    dots_list <- list(...)
    dots_list <- dots_list[!vapply(dots_list, is.null, logical(1))]

    single_network_args <- c(network_params, dots_list)

    network_table <- thisutils::parallelize_fun(
      x = targets,
      fun = function(x) {
        thisutils::invoke_fun(
          single_network,
          c(
            list(
              matrix = object,
              regulators = regulators,
              target = x
            ),
            single_network_args
          )
        )
      },
      cores = validated_params$cores,
      verbose = verbose
    ) |>
      purrr::list_rbind() |>
      network_format(abs_weight = FALSE)

    thisutils::log_message(
      "Inferring network done",
      message_type = "success",
      verbose = verbose
    )
    network_info <- data.frame(
      Edges = nrow(network_table),
      Regulators = length(unique(network_table$regulator)),
      Targets = length(unique(network_table$target))
    )
    thisutils::log_message(
      "Network information:\n",
      network_info,
      text_color = "grey",
      timestamp_style = FALSE,
      verbose = verbose
    )

    return(network_table)
  }
)
