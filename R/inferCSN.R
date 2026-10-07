#' @title inferring cell-type specific gene regulatory network
#'
#' @md
#' @param object Numeric expression matrix with cells in rows and genes in
#' columns.
#' @param pseudotime Optional pseudotime vector or branch matrix for either method.
#'   Regulators in earlier states predict targets in later states; tied states
#'   are averaged and shared branch transitions are counted once.
#' @param regulators,targets Optional gene subsets.
#' @param max_support_size Optional support-size limit.
#' @param lag_fraction Fractional lag used when `lag_steps` is `NULL`.
#' @param lag_steps Optional integer lag. L0Learn centers and scales lagged
#'   expression within each branch and defaults to fitting without an intercept.
#' @param cores Number of inference workers.
#' @param verbose Whether to report progress.
#' @param method `greedy_l0` (default), or L0Learn with the `L0`, `L0L1`,
#'   or `L0L2` penalty.
#' @param ... Arguments passed to the method.
#'
#' @return A data frame containing exactly `regulator`, `target`, and `weight`.
#'
#' @docType methods
#' @rdname inferCSN
#' @export
methods::setGeneric(
  name = "inferCSN",
  signature = "object",
  def = function(
    object,
    pseudotime = NULL,
    regulators = NULL,
    targets = NULL,
    max_support_size = NULL,
    lag_fraction = 0.05,
    lag_steps = NULL,
    cores = 1,
    verbose = TRUE,
    method = c("greedy_l0", "L0", "L0L1", "L0L2"),
    ...
  ) {
    standardGeneric("inferCSN")
  }
)

infercsn_method <- function(
  object,
  pseudotime = NULL,
  regulators = NULL,
  targets = NULL,
  max_support_size = NULL,
  lag_fraction = 0.05,
  lag_steps = NULL,
  cores = 1,
  verbose = TRUE,
  method = c("greedy_l0", "L0", "L0L1", "L0L2"),
  ...
) {
  method <- match.arg(method)
  dots <- list(...)
  if (identical(method, "greedy_l0") && length(dots)) {
    thisutils::log_message(
      sprintf(
        "Unused matrix-inference argument%s: %s",
        if (length(dots) == 1L) "" else "s",
        paste(names(dots) %|||% rep("<unnamed>", length(dots)), collapse = ", ")
      ),
      message_type = "error"
    )
  }
  thisutils::log_message(
    "Inferring network for {.cls {class(object)}}...",
    verbose = verbose
  )
  validated <- validate_parameters(
    matrix = object,
    pseudotime = pseudotime,
    regulators = regulators,
    targets = targets,
    max_support_size = max_support_size,
    lag_fraction = lag_fraction,
    lag_steps = lag_steps,
    cores = cores,
    method = method,
    dots = dots,
    verbose = verbose
  )
  pseudotime <- validated$pseudotime
  if (!identical(method, "greedy_l0")) {
    regulators <- validated$regulators %|||% colnames(object)
    targets <- validated$targets %|||% colnames(object)
    if (length(regulators) < 2L) {
      thisutils::log_message("L0Learn network inference requires at least two regulators.",
        message_type = "error"
      )
    }
    names(targets) <- targets
    parameters <- validated
    parameters$cores <- 1L
    fits <- thisutils::parallelize_fun(
      x = targets,
      fun = function(target) {
        fit_l0learn_target(
          object, regulators, target, parameters, method, verbose, ...
        )
      },
      cores = validated$cores, verbose = verbose, clean_result = TRUE
    )
    network_table <- do.call(rbind, fits)
    if (is.null(network_table)) {
      network_table <- data.frame(regulator = character(), target = character(), weight = numeric())
    } else {
      network_table <- network_format(network_table, abs_weight = FALSE)
    }
  } else {
    network_table <- infer_network(
      expression = as.matrix(object),
      pseudotime = pseudotime,
      gene_names = as.character(colnames(object)),
      params = list(
        min_improvement = 1e-10,
        pseudotime_lag_fraction = validated$lag_fraction,
        pseudotime_lag_steps = validated$lag_steps,
        regulators = validated$regulators,
        targets = validated$targets,
        max_support_size = validated$max_support_size,
        cores = validated$cores
      )
    )
    network_table <- network_table[
      is.finite(network_table$weight) & network_table$weight != 0,
      c("regulator", "target", "weight"),
      drop = FALSE
    ]
  }

  thisutils::log_message(
    "Inferring network done",
    message_type = "success",
    verbose = verbose
  )
  thisutils::log_message(
    "Network information:\n",
    data.frame(
      Edges = nrow(network_table),
      Regulators = length(unique(network_table$regulator)),
      Targets = length(unique(network_table$target))
    ),
    text_color = "grey",
    timestamp_style = FALSE,
    verbose = verbose
  )
  network_table
}

#' @rdname inferCSN
#' @export
#' @examples
#' data(example_matrix)
#' data(example_meta_data)
#' network_table <- inferCSN(
#'   example_matrix,
#'   pseudotime = example_meta_data$pseudotime
#' )
#' head(network_table)
#'
#' inferCSN(
#'   example_matrix,
#'   regulators = c("g1", "g2"),
#'   targets = c("g3", "g4")
#' )
methods::setMethod(
  f = "inferCSN",
  signature = methods::signature(object = "matrix"),
  definition = infercsn_method
)

#' @rdname inferCSN
#' @export
methods::setMethod(
  f = "inferCSN",
  signature = methods::signature(object = "sparseMatrix"),
  definition = infercsn_method
)
