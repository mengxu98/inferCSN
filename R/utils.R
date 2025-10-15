.check_parameters <- function(
    matrix,
    penalty,
    cross_validation,
    seed,
    n_folds,
    r_squared_threshold,
    regulators,
    targets,
    verbose,
    cores,
    ...) {
  thisutils::log_message(
    "Checking input parameters...",
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

  match.arg(penalty, c("L0", "L0L1", "L0L2"))

  if (!is.numeric(seed)) {
    seed <- 1
    thisutils::log_message(
      "initialize random seed to 1",
      message_type = "warning",
      verbose = verbose
    )
  }

  if (r_squared_threshold < 0 || r_squared_threshold > 1) {
    thisutils::log_message(
      "Please set {.arg r_squared_threshold} between: [0, 1]",
      message_type = "error"
    )
  }

  if (!is.null(regulators)) {
    intersect_regulators <- intersect(regulators, colnames(matrix))
    if (length(intersect_regulators) < 2) {
      thisutils::log_message(
        "The input genes must contain at least 2 regulators",
        message_type = "error"
      )
    }

    if (length(intersect_regulators) < length(regulators)) {
      thisutils::log_message(
        "{.val {length(intersect_regulators)}} out of {.val {length(regulators)}} candidate regulator{?s} are in the input matrix",
        message_type = "warning",
        verbose = verbose
      )
    } else {
      thisutils::log_message(
        "Using {.val {length(intersect_regulators)}} regulator{?s}",
        verbose = verbose
      )
    }
  }

  if (!is.null(targets)) {
    intersect_targets <- intersect(targets, colnames(matrix))
    if (length(intersect_targets) == 0) {
      thisutils::log_message(
        "The input genes must contain at least 1 target",
        message_type = "error"
      )
    }

    if (length(intersect_targets) < length(targets)) {
      thisutils::log_message(
        "{.val {length(intersect_targets)}} out of {.val {length(targets)}} candidate target{?s} are in the input matrix",
        message_type = "warning",
        verbose = verbose
      )
    } else {
      thisutils::log_message(
        "Using {.val {length(intersect_targets)}} target{?s}",
        verbose = verbose
      )
    }
  }

  if (!is.numeric(cores) || cores < 1) {
    thisutils::log_message(
      "{.arg cores} should be a positive integer, initialize it to 1",
      message_type = "warning",
      verbose = verbose
    )
  }

  thisutils::log_message(
    "Using {.arg {penalty}} sparse regression model",
    verbose = verbose
  )
  if (cross_validation) {
    thisutils::log_message(
      "Using cross validation, and setting {.val {n_folds}} fold{?s}",
      verbose = verbose
    )
  }
}
