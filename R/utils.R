validate_parameters <- function(
  matrix,
  pseudotime,
  regulators,
  targets,
  max_support_size,
  lag_fraction,
  lag_steps,
  cores,
  method,
  dots,
  verbose
) {
  controls <- NULL
  if (!identical(method, "greedy_l0")) {
    if ("penalty" %in% names(dots)) {
      thisutils::log_message("Select the L0Learn penalty with `method`, not `penalty`.", message_type = "error")
    }
    if (!requireNamespace("L0Learn", quietly = TRUE)) {
      thisutils::log_message("Install L0Learn to use method = 'L0', 'L0L1', or 'L0L2'.", message_type = "error")
    }
    controls <- list(
      cross_validation = dots$cross_validation %|||% FALSE,
      seed = dots$seed %|||% 1,
      n_folds = dots$n_folds %|||% 5,
      r_squared_threshold = dots$r_squared_threshold %|||% 0
    )
    if (!is.logical(controls$cross_validation) || length(controls$cross_validation) != 1L ||
      is.na(controls$cross_validation)) {
      thisutils::log_message("`cross_validation` must be TRUE or FALSE.", message_type = "error")
    }
    if (!is.numeric(controls$seed) || length(controls$seed) != 1L || !is.finite(controls$seed)) {
      thisutils::log_message("`seed` must be one finite number.", message_type = "error")
    }
    if (!is.numeric(controls$n_folds) || length(controls$n_folds) != 1L || !is.finite(controls$n_folds) ||
      controls$n_folds < 2 || controls$n_folds != as.integer(controls$n_folds)) {
      thisutils::log_message("`n_folds` must be an integer >= 2.", message_type = "error")
    }
    if (!is.numeric(controls$r_squared_threshold) || length(controls$r_squared_threshold) != 1L ||
      !is.finite(controls$r_squared_threshold) || controls$r_squared_threshold < 0 || controls$r_squared_threshold > 1) {
      thisutils::log_message("`r_squared_threshold` must be between 0 and 1.", message_type = "error")
    }
  }

  if (length(dim(matrix)) != 2L) {
    thisutils::log_message("`object` must be a two-dimensional matrix.", message_type = "error")
  }
  if (is.null(colnames(matrix)) || anyNA(colnames(matrix)) ||
    any(!nzchar(colnames(matrix))) || anyDuplicated(colnames(matrix))) {
    thisutils::log_message("`object` must have unique, non-empty gene names as column names.", message_type = "error")
  }

  if (!is.numeric(cores) || length(cores) != 1L || !is.finite(cores) ||
    cores < 1 || cores != as.integer(cores)) {
    thisutils::log_message("`cores` must be one positive integer.", message_type = "error")
  }
  if (!is.numeric(lag_fraction) || length(lag_fraction) != 1L ||
    !is.finite(lag_fraction) || lag_fraction <= 0 || lag_fraction > 1) {
    thisutils::log_message("`lag_fraction` must be one number in (0, 1].", message_type = "error")
  }
  if (is.null(lag_steps)) {
    lag_steps <- 0L
  } else if (!is.numeric(lag_steps) || length(lag_steps) != 1L ||
    !is.finite(lag_steps) || lag_steps < 1 ||
    lag_steps != as.integer(lag_steps)) {
    thisutils::log_message("`lag_steps` must be `NULL` or one positive integer.", message_type = "error")
  }

  max_support_size <- validate_max_support_size(max_support_size)

  validate_genes <- function(requested, label) {
    if (is.null(requested)) {
      return(NULL)
    }
    requested <- unique(as.character(requested))
    requested <- requested[!is.na(requested) & nzchar(requested)]
    present <- intersect(requested, colnames(matrix))
    missing <- setdiff(requested, colnames(matrix))
    if (!length(present)) {
      thisutils::log_message(sprintf("None of the requested %s are present in `object`.", label), message_type = "error")
    }
    if (length(missing)) {
      thisutils::log_message(
        sprintf(
          "Ignoring %d requested %s absent from `object`: %s",
          length(missing),
          label,
          paste(missing, collapse = ", ")
        ),
        message_type = "warning", verbose = verbose
      )
    }
    present
  }

  list(
    pseudotime = validate_pseudotime(pseudotime, nrow(matrix)),
    regulators = validate_genes(regulators, "regulators"),
    targets = validate_genes(targets, "targets"),
    max_support_size = as.integer(max_support_size),
    lag_fraction = as.numeric(lag_fraction),
    lag_steps = as.integer(lag_steps),
    cores = as.integer(cores),
    controls = controls
  )
}

`%|||%` <- thisutils::`%|||%`

validate_pseudotime <- function(pseudotime, n_cells) {
  if (is.null(pseudotime)) {
    return(matrix(numeric(0L), nrow = n_cells, ncol = 0L))
  }
  pseudotime <- if (is.data.frame(pseudotime) || is.matrix(pseudotime)) {
    as.matrix(pseudotime)
  } else {
    matrix(pseudotime, ncol = 1L)
  }
  if (!is.numeric(pseudotime)) {
    thisutils::log_message("`pseudotime` must be numeric.", message_type = "error")
  }
  if (nrow(pseudotime) != n_cells) {
    thisutils::log_message(
      "`pseudotime` must contain one row (or one vector value) per cell.",
      message_type = "error"
    )
  }
  storage.mode(pseudotime) <- "double"
  pseudotime
}
