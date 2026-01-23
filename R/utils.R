check_parameters <- function(
    matrix,
    penalty,
    cross_validation,
    seed,
    n_folds,
    subsampling_method,
    subsampling_ratio,
    r_squared_threshold,
    regulators,
    targets,
    cores,
    verbose,
    ...) {
  thisutils::log_message(
    "Checking parameters...",
    message_type = "running",
    verbose = verbose
  )

  if (length(dim(matrix)) != 2) {
    thisutils::log_message(
      "The input matrix must be a two-dimensional matrix",
      message_type = "error"
    )
  }

  if (is.null(colnames(matrix))) {
    thisutils::log_message(
      "The input matrix must contain the names of the genes as colnames",
      message_type = "error"
    )
  }

  param_names <- c(
    "penalty", "subsampling_method", "subsampling_ratio",
    "n_folds", "cores", "seed", "r_squared_threshold", "cross_validation"
  )
  params <- mget(param_names, ifnotfound = list(NULL))

  validated_params <- list()

  validated_params$penalty <- match.arg(params$penalty, choices = c("L0", "L0L1", "L0L2"))

  validated_params$subsampling_method <- match.arg(
    params$subsampling_method,
    choices = c("sample", "meta_cells", "pseudobulk")
  )

  if (!is.numeric(params$subsampling_ratio) || params$subsampling_ratio <= 0 || params$subsampling_ratio > 1) {
    thisutils::log_message(
      "Please set {.arg subsampling_ratio} between: (0, 1]",
      message_type = "error"
    )
  }
  validated_params$subsampling_ratio <- params$subsampling_ratio

  if (!is.numeric(params$n_folds) || params$n_folds < 2 || params$n_folds != as.integer(params$n_folds)) {
    thisutils::log_message(
      "Please set {.arg n_folds} as an integer >= 2",
      message_type = "error"
    )
  }
  validated_params$n_folds <- as.integer(params$n_folds)

  if (!is.numeric(params$cores) || params$cores < 1 || params$cores != as.integer(params$cores)) {
    thisutils::log_message(
      "Please set {.arg cores} as an integer >= 1",
      message_type = "error"
    )
  }
  validated_params$cores <- as.integer(params$cores)

  if (!is.numeric(params$seed)) {
    thisutils::log_message(
      "Initialize random seed to 1",
      message_type = "warning",
      verbose = verbose
    )
    validated_params$seed <- 1
  } else {
    validated_params$seed <- params$seed
  }

  if (params$r_squared_threshold < 0 || params$r_squared_threshold > 1) {
    thisutils::log_message(
      "Please set {.arg r_squared_threshold} between: [0, 1]",
      message_type = "error"
    )
  }
  validated_params$r_squared_threshold <- params$r_squared_threshold

  if (!is.logical(params$cross_validation)) {
    thisutils::log_message(
      "Please set {.arg cross_validation} as a logical value",
      message_type = "error"
    )
  }
  validated_params$cross_validation <- params$cross_validation

  thisutils::log_message(
    "Using {.pkg {validated_params$penalty}} sparse regression model",
    verbose = verbose
  )

  validated_params$regulators <- NULL
  if (!is.null(regulators)) {
    intersect_regulators <- intersect(regulators, colnames(matrix))
    missing_regulators <- setdiff(regulators, colnames(matrix))

    if (length(intersect_regulators) < 2) {
      thisutils::log_message(
        "The input genes must contain at least 2 regulators",
        message_type = "error"
      )
    }

    if (length(missing_regulators) > 0) {
      thisutils::log_message(
        "The following regulator{?s} {?is/are} not in the matrix: {.val {missing_regulators}}",
        message_type = "warning",
        verbose = verbose
      )
    }

    if (length(intersect_regulators) < length(regulators)) {
      thisutils::log_message(
        "{.val {length(intersect_regulators)}} out of {.val {length(regulators)}} regulator{?s} are in the input matrix",
        message_type = "warning",
        verbose = verbose
      )
    } else {
      thisutils::log_message(
        "Using {.val {length(intersect_regulators)}} regulator{?s}",
        verbose = verbose
      )
    }

    validated_params$regulators <- intersect_regulators
  }

  validated_params$targets <- NULL
  if (!is.null(targets)) {
    intersect_targets <- intersect(targets, colnames(matrix))
    missing_targets <- setdiff(targets, colnames(matrix))

    if (length(intersect_targets) == 0) {
      thisutils::log_message(
        "The input genes must contain at least 1 target",
        message_type = "error"
      )
    }

    if (length(missing_targets) > 0) {
      thisutils::log_message(
        "The following target{?s} {?is/are} not in the matrix: {.val {missing_targets}}",
        message_type = "warning",
        verbose = verbose
      )
    }

    if (length(intersect_targets) < length(targets)) {
      thisutils::log_message(
        "{.val {length(intersect_targets)}} out of {.val {length(targets)}} target{?s} are in the input matrix",
        message_type = "warning",
        verbose = verbose
      )
    } else {
      thisutils::log_message(
        "Using {.val {length(intersect_targets)}} target{?s}",
        verbose = verbose
      )
    }

    validated_params$targets <- intersect_targets
  }

  if (validated_params$cross_validation) {
    thisutils::log_message(
      "Using {.val {validated_params$n_folds}} fold{?s} cross validation",
      verbose = verbose
    )
  }

  return(validated_params)
}
