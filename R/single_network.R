#' @title Construct network for single target gene
#'
#' @inheritParams inferCSN
#' @param matrix An expression matrix.
#' @param regulators Candidate regulator genes.
#' @param target The target gene.
#'
#' @param pseudotime Optional pseudotime vector or branch matrix passed to
#' [inferCSN()].
#' @param max_support_size Optional support-size cap passed to [inferCSN()].
#' @param lag_fraction Fractional state lag passed to [inferCSN()].
#' @param lag_steps Optional integer state lag passed to [inferCSN()].
#' @param cores Number of inference workers.
#'
#' @return A data frame containing only selected edges for the requested target.
#' The data frame has three columns: regulator, target, and weight.
#'
#' @export
#' @examples
#' data(example_matrix)
#' head(
#'   single_network(
#'     example_matrix,
#'     regulators = colnames(example_matrix),
#'     target = "g1"
#'   )
#' )
#' single_network(
#'   example_matrix,
#'   regulators = c("g1", "g2", "g3"),
#'   target = "g1"
#' )
single_network <- function(
  matrix,
  regulators,
  target,
  pseudotime = NULL,
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
  if (!identical(method, "greedy_l0")) {
    parameters <- validate_parameters(
      matrix, pseudotime, regulators, target,
      max_support_size, lag_fraction, lag_steps, cores, method, dots, verbose
    )
    return(fit_l0learn_target(
      matrix, parameters$regulators, target,
      parameters, method, verbose, ...
    ))
  }
  if (length(dots)) {
    thisutils::log_message("L0Learn fitting arguments require method = 'L0', 'L0L1', or 'L0L2'.", message_type = "error")
  }
  regulators <- setdiff(regulators, target)
  if (length(regulators) < 1) {
    thisutils::log_message(
      "No candidate regulators found when modeling: {.val {target}}",
      message_type = "warning",
      verbose = verbose
    )
    return(data.frame(
      regulator = character(),
      target = character(),
      weight = numeric()
    ))
  }
  inferCSN(
    object = matrix,
    pseudotime = pseudotime,
    regulators = regulators,
    targets = target,
    max_support_size = max_support_size,
    lag_fraction = lag_fraction,
    lag_steps = lag_steps,
    cores = cores,
    verbose = verbose
  )
}

fit_l0learn_target <- function(
  matrix, regulators, target, parameters, method, verbose, ...
) {
  controls <- parameters$controls
  dots <- list(...)
  dots[c("cross_validation", "seed", "n_folds", "r_squared_threshold")] <- NULL
  regulators <- setdiff(regulators, target)
  if (length(regulators) < 2L) {
    thisutils::log_message(
      "Less than 2 regulators found when modeling: {.val {target}}",
      message_type = "warning", verbose = verbose
    )
    return(data.frame(regulator = character(), target = character(), weight = numeric()))
  }
  if (ncol(parameters$pseudotime)) {
    genes <- c(regulators, target)
    aligned <- prepare_lagged_expression(
      as.matrix(matrix[, genes, drop = FALSE]),
      parameters$pseudotime, parameters$lag_fraction, parameters$lag_steps, parameters$cores
    )
    x <- aligned$x[, regulators, drop = FALSE]
    y <- aligned$y[, target]
    if (aligned$lagged) dots$intercept <- dots$intercept %|||% FALSE
  } else {
    x <- matrix[, regulators]
    y <- matrix[, target]
  }
  fit_args <- c(list(
    x = x, y = y, penalty = method,
    maxSuppSize = if (parameters$max_support_size == 0L) {
      ncol(x)
    } else {
      min(parameters$max_support_size, ncol(x))
    }
  ), dots)
  log_warning <- function(w) {
    thisutils::log_message(conditionMessage(w), message_type = "warning", verbose = verbose)
    invokeRestart("muffleWarning")
  }
  fit <- tryCatch(withCallingHandlers(
    if (controls$cross_validation) {
      do.call(
        L0Learn::L0Learn.cvfit,
        c(fit_args, list(nFolds = as.integer(controls$n_folds), seed = controls$seed))
      )
    } else {
      do.call(L0Learn::L0Learn.fit, fit_args)
    },
    warning = log_warning
  ), error = identity)
  if (controls$cross_validation && inherits(fit, "error")) {
    thisutils::log_message(
      "Cross-validation error, setting {.arg cross_validation} to {.pkg FALSE} and re-train",
      message_type = "warning", verbose = verbose
    )
    fit <- tryCatch(withCallingHandlers(do.call(L0Learn::L0Learn.fit, fit_args),
      warning = log_warning
    ), error = identity)
  }
  coefficients <- rep(0, ncol(x))
  if (!inherits(fit, "error")) {
    if (inherits(fit, "L0LearnCV")) {
      gamma_index <- which.min(sapply(fit$cvMeans, min))
      gamma <- fit$fit$gamma[gamma_index]
      lambda_index <- which.min(fit$cvMeans[[gamma_index]])
      lambda <- fit$fit$lambda[[gamma_index]][lambda_index]
    } else {
      fit_info <- data.frame(
        lambda = unlist(fit$lambda),
        gamma = rep(fit$gamma, lengths(fit$lambda)), suppSize = unlist(fit$suppSize)
      )
      index <- which.max(fit_info$suppSize)
      lambda <- fit_info$lambda[index]
      gamma <- fit_info$gamma[index]
    }
    pred_y <- as.numeric(predict(fit, newx = x, lambda = lambda, gamma = gamma))
    if (thisutils::r_square(y, pred_y) > controls$r_squared_threshold) {
      coefficients <- as.vector(coef(fit, lambda = lambda, gamma = gamma))
      fitted_model <- if (inherits(fit, "L0LearnCV")) fit$fit else fit
      if (fitted_model$settings$intercept) coefficients <- coefficients[-1]
      coefficients <- thisutils::normalization(coefficients, method = "unit_vector", ...)
      if (length(coefficients) != ncol(x)) {
        coefficients <- rep(0, ncol(x))
      }
    }
  }
  if (inherits(fit, "error")) {
    thisutils::log_message("Fitting failed for {.val {target}}: {conditionMessage(fit)}",
      message_type = "warning", verbose = verbose
    )
  }
  network_format(data.frame(regulator = regulators, target = target, weight = coefficients), abs_weight = FALSE)
}
