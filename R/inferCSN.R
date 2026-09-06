#' @title inferring cell-type specific gene regulatory network
#'
#' @md
#' @description Fits greedy-L0 models for static or pseudotime-ordered expression data.
#'
#' @param object Numeric expression matrix with cells in rows and genes in
#' columns.
#' @param pseudotime Optional pseudotime vector or branch matrix.
#' @param regulators,targets Optional gene subsets.
#' @param max_support_size Optional support-size limit.
#' @param lag_fraction Fractional lag used when `lag_steps` is `NULL`.
#' @param lag_steps Optional integer lag.
#' @param cores Number of inference workers.
#' @param verbose Whether to report progress.
#' @param ... Additional method arguments.
#'
#' @return A data frame containing exactly `regulator`, `target`, and `weight`.
#' @details Signed ordinal weights group descending deletion evidence against
#' each group's maximum within 1e-12 * (1 + abs(maximum)). This fixed numerical
#' rule leaves support, fitted coefficients and raw deletion evidence unchanged.
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
  ...
) {
  dots <- list(...)
  if (length(dots)) {
    stop(
      sprintf(
        "Unused matrix-inference argument%s: %s",
        if (length(dots) == 1L) "" else "s",
        paste(names(dots) %|||% rep("<unnamed>", length(dots)), collapse = ", ")
      ),
      call. = FALSE
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
    verbose = verbose
  )

  gene_names <- colnames(object)
  expression <- t(as.matrix(object))
  pseudotime <- validated$pseudotime
  n_cells <- ncol(expression)
  if (is.null(pseudotime)) {
    pseudotime <- matrix(numeric(0L), nrow = n_cells, ncol = 0L)
  } else if (is.data.frame(pseudotime) || is.matrix(pseudotime)) {
    pseudotime <- as.matrix(pseudotime)
  } else {
    pseudotime <- matrix(pseudotime, ncol = 1L)
  }
  if (!is.numeric(pseudotime)) {
    stop("`pseudotime` must be numeric.", call. = FALSE)
  }
  if (nrow(pseudotime) != n_cells) {
    stop(
      "`pseudotime` must contain one row (or one vector value) per cell.",
      call. = FALSE
    )
  }
  storage.mode(pseudotime) <- "double"
  network_table <- infer_network(
    expression = expression,
    pseudotime = pseudotime,
    gene_names = as.character(gene_names),
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
