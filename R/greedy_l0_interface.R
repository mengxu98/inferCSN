run_network_inference <- function(expression, pseudotime, gene_names, params) {
  if (inherits(expression, "sparseMatrix")) {
    expression <- as.matrix(expression)
    expression <- t(expression)
  }
  if (!is.matrix(expression)) {
    expression <- as.matrix(expression)
  }
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
  infer_network(
    expression, as.character(gene_names), pseudotime, params
  )
}
