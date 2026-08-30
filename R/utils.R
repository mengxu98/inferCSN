validate_infercsn_parameters <- function(
  matrix,
  pseudotime,
  regulators,
  targets,
  max_support_size,
  lag_fraction,
  lag_steps,
  cores,
  verbose
) {
  thisutils::log_message(
    "Checking parameters...",
    message_type = "running",
    verbose = verbose
  )

  if (length(dim(matrix)) != 2L) {
    stop("`object` must be a two-dimensional matrix.", call. = FALSE)
  }
  if (is.null(colnames(matrix)) || anyNA(colnames(matrix)) ||
    any(!nzchar(colnames(matrix))) || anyDuplicated(colnames(matrix))) {
    stop("`object` must have unique, non-empty gene names as column names.", call. = FALSE)
  }

  if (!is.numeric(cores) || length(cores) != 1L || !is.finite(cores) ||
    cores < 1 || cores != as.integer(cores)) {
    stop("`cores` must be one positive integer.", call. = FALSE)
  }
  if (!is.numeric(lag_fraction) || length(lag_fraction) != 1L ||
    !is.finite(lag_fraction) || lag_fraction <= 0 || lag_fraction > 1) {
    stop("`lag_fraction` must be one number in (0, 1].", call. = FALSE)
  }
  if (is.null(lag_steps)) {
    lag_steps <- 0L
  } else if (!is.numeric(lag_steps) || length(lag_steps) != 1L ||
    !is.finite(lag_steps) || lag_steps < 1 ||
    lag_steps != as.integer(lag_steps)) {
    stop("`lag_steps` must be `NULL` or one positive integer.", call. = FALSE)
  }

  max_support_size <- validate_max_support_size(max_support_size)
  if (is.null(max_support_size)) {
    max_support_size <- 0L
  }

  validate_genes <- function(requested, label) {
    if (is.null(requested)) {
      return(NULL)
    }
    requested <- unique(as.character(requested))
    requested <- requested[!is.na(requested) & nzchar(requested)]
    present <- intersect(requested, colnames(matrix))
    missing <- setdiff(requested, colnames(matrix))
    if (!length(present)) {
      stop(sprintf("None of the requested %s are present in `object`.", label), call. = FALSE)
    }
    if (length(missing)) {
      warning(
        sprintf(
          "Ignoring %d requested %s absent from `object`: %s",
          length(missing),
          label,
          paste(missing, collapse = ", ")
        ),
        call. = FALSE
      )
    }
    present
  }

  list(
    pseudotime = pseudotime,
    regulators = validate_genes(regulators, "regulators"),
    targets = validate_genes(targets, "targets"),
    max_support_size = as.integer(max_support_size),
    lag_fraction = as.numeric(lag_fraction),
    lag_steps = as.integer(lag_steps),
    cores = as.integer(cores)
  )
}

`%||%` <- function(x, y) {
  if (is.null(x) || !length(x)) y else x
}
