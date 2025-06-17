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
    "Checking input parameters.",
    verbose = verbose
  )

  if (length(dim(matrix)) != 2) {
    thisutils::log_message(
      "The input matrix must be a two-dimensional matrix.",
      message_type = "error"
    )
  }

  if (is.null(colnames(matrix))) {
    thisutils::log_message(
      "The input matrix must contain the names of the genes as colnames.",
      message_type = "error"
    )
  }

  match.arg(penalty, c("L0", "L0L1", "L0L2"))

  if (!is.numeric(seed)) {
    seed <- 1
    thisutils::log_message(
      "initialize random seed to 1.",
      message_type = "warning",
      verbose = verbose
    )
  }

  if (r_squared_threshold < 0 || r_squared_threshold > 1) {
    thisutils::log_message(
      "Please set 'r_squared_threshold' value between: [0, 1].",
      message_type = "error"
    )
  }

  if (!is.null(regulators)) {
    intersect_regulators <- intersect(regulators, colnames(matrix))
    if (length(intersect_regulators) < 2) {
      thisutils::log_message(
        "The input genes must contain at least 2 regulator.",
        message_type = "error"
      )
    }

    if (length(intersect_regulators) < length(regulators)) {
      thisutils::log_message(
        length(intersect_regulators), " out of ",
        length(regulators), " candidate regulators are in the input matrix.",
        message_type = "warning",
        verbose = verbose
      )
    } else {
      thisutils::log_message(
        "Using ", length(intersect_regulators), " regulator(s).",
        verbose = verbose
      )
    }
  }

  if (!is.null(targets)) {
    intersect_targets <- intersect(targets, colnames(matrix))
    if (length(intersect_targets) == 0) {
      thisutils::log_message(
        "The input genes must contain at least 1 target.",
        message_type = "error"
      )
    }

    if (length(intersect_targets) < length(targets)) {
      thisutils::log_message(
        length(intersect_targets), " out of ",
        length(targets), " candidate targets are in the input matrix.",
        message_type = "warning",
        verbose = verbose
      )
    } else {
      thisutils::log_message(
        "Using ", length(intersect_targets), " target(s).",
        verbose = verbose
      )
    }
  }

  if (!is.numeric(cores) || cores < 1) {
    thisutils::log_message(
      "`cores` should be a positive integer, initialize it to 1.",
      message_type = "warning",
      verbose = verbose
    )
  }

  thisutils::log_message(
    "Using ", penalty, " sparse regression model.",
    verbose = verbose
  )
  if (cross_validation) {
    thisutils::log_message(
      "Using cross validation, and setting ", n_folds, " folds.",
      verbose = verbose
    )
  }
}

#' @title Extracts a specific solution in the regularization path
#'
#' @inheritParams fit_srm
#' @param object The output of \code{\link{sparse_regression}}.
#' @param lambda The value of lambda at which to extract the solution.
#' @param gamma The value of gamma at which to extract the solution.
#' @param ... Other parameters
#'
#' @method coef srm
#'
#' @return Return the specific solution
#' @export
coef.srm <- function(
    object,
    lambda = NULL,
    gamma = NULL,
    regulators_num = NULL,
    ...) {
  if (!is.null(regulators_num) && !is.null(lambda)) {
    stop("If 'regulators_num' is provided to 'coef' only 'gamma' can also be provided.")
  }

  if (is.null(lambda) && is.null(gamma) && is.null(regulators_num)) {
    t <- do.call(cbind, object$beta)
    if (object$settings$intercept) {
      intercepts <- unlist(object$a0)
      t <- rbind(intercepts, t)
    }
    return(t)
  }

  if (is.null(gamma)) {
    gamma <- object$gamma[1]
  }

  diff_gamma <- abs(object$gamma - gamma)
  gamma_index <- which(diff_gamma == min(diff_gamma))

  indices <- NULL
  if (!is.null(lambda)) {
    diff_lambda <- abs(lambda - object$lambda[[gamma_index]])
    indices <- which(diff_lambda == min(diff_lambda))
  } else if (!is.null(regulators_num)) {
    diff_regulators_num <- abs(
      regulators_num - object$suppSize[[gamma_index]]
    )
    indices <- which(
      diff_regulators_num == min(diff_regulators_num)
    )
  } else {
    indices <- seq_along(object$lambda[[gamma_index]])
  }

  if (object$settings$intercept) {
    t <- rbind(
      object$a0[[gamma_index]][indices],
      object$beta[[gamma_index]][, indices, drop = FALSE]
    )
    rownames(t) <- c(
      "Intercept",
      paste0(rep("V", object$p), 1:object$p)
    )
  } else {
    t <- object$beta[[gamma_index]][, indices, drop = FALSE]
    rownames(t) <- paste0(rep("V", object$p), 1:object$p)
  }

  return(t)
}

#' @rdname coef.srm
#' @method coef srm_cv
#' @export
coef.srm_cv <- function(
    object,
    lambda = NULL,
    gamma = NULL,
    ...) {
  return(
    coef.srm(
      object$fit, lambda, gamma, ...
    )
  )
}

#' @title Prints a summary of `sparse_regression`
#'
#' @param x The output of \code{\link{sparse_regression}}.
#' @param ... Other parameters
#'
#' @method print srm
#'
#' @md
#' @return Return information of `sparse_regression`
#' @export
print.srm <- function(x, ...) {
  gammas <- rep(x$gamma, times = lapply(x$lambda, length))
  return(
    data.frame(
      lambda = unlist(x["lambda"]),
      gamma = gammas,
      suppSize = unlist(x["suppSize"]),
      row.names = NULL
    )
  )
}

#' @rdname print.srm
#' @method print srm_cv
#' @export
print.srm_cv <- function(x, ...) {
  return(
    print.srm(x$fit)
  )
}

#' @title Predicts response for a given sample
#'
#' @param object The output of sparse_regression.
#' @param newx A matrix on which predictions are made. The matrix should have p columns
#' @param lambda The value of lambda to use for prediction.
#' A summary of the lambdas in the regularization path can be obtained using \code{\link{print.srm}}.
#' @param gamma The value of gamma to use for prediction.
#' A summary of the gammas in the regularization path can be obtained using \code{\link{print.srm}}.
#' @param ... Other parameters
#'
#' @method predict srm
#'
#' @details
#' If both lambda and gamma are not supplied,
#' then a matrix of predictions for all the solutions in the regularization path is returned.
#' If lambda is supplied but gamma is not, the smallest value of gamma is used.
#' In case of logistic regression, probability values are returned.
#'
#' @return Return predict value
#' @export
predict.srm <- function(
    object,
    newx,
    lambda = NULL,
    gamma = NULL,
    ...) {
  beta <- coef.srm(object, lambda, gamma)
  if (object$settings$intercept) {
    # add a column of ones for the intercept
    x <- cbind(1, newx)
  } else {
    x <- newx
  }
  prediction <- x %*% beta
  if (object$loss == "Logistic") {
    prediction <- 1 / (1 + exp(-prediction))
  }

  return(prediction)
}

#' @rdname predict.srm
#' @method predict srm_cv
#' @export
predict.srm_cv <- function(
    object,
    newx,
    lambda = NULL,
    gamma = NULL,
    ...) {
  return(
    predict.srm(
      object$fit, newx, lambda, gamma, ...
    )
  )
}

.is_scalar <- function(x) {
  is.atomic(x) && length(x) == 1L && !is.character(x) && Im(x) == 0 && !is.nan(x) && !is.na(x)
}
