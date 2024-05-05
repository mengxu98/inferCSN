#' @title Check input parameters
#'
#' @param matrix An expression matrix, cells by genes
#' @inheritParams inferCSN
#'
#' @return No return value, called for check input parameters
#' @export
check.parameters <- function(
    matrix,
    penalty,
    algorithm,
    cross_validation,
    seed,
    n_folds,
    k_folds,
    r_threshold,
    regulators,
    targets,
    regulators_num,
    verbose,
    cores) {
  if (verbose) {
    message("Checking input parameters.")
  }

  if (!is.matrix(matrix) && !is.array(matrix) && (length(dim(matrix)) != 2)) {
    stop(
      "`matrix` must be a two-dimensional matrix,
      where each column corresponds to a gene and each row corresponds to a sample/cell."
    )
  }

  if (is.null(colnames(matrix))) {
    stop("`matrix` must contain the names of the genes as colnames.")
  }

  # Check the penalty term of the regression model
  match.arg(penalty, c("L0", "L0L2"))

  # Check the algorithm of the regression model
  match.arg(algorithm, c("CD", "CDPSI"))

  if (!is.numeric(seed)) {
    seed <- 1
    if (verbose) {
      warning("`seed` is not a valid integer, initialize to 1.")
    }
  }

  if (!is.null(k_folds)) {
    if (!(is.numeric(k_folds) && k_folds > 0 && k_folds <= 1)) {
      stop("Please set `k_folds` value between: (0, 1].")
    }
  }

  if (!is.null(targets)) {
    targets_intersect <- intersect(targets, colnames(matrix))
    if (length(targets_intersect) == 0) {
      stop("The input genes must contain at least 1 target.")
    }

    if (length(targets_intersect) < length(targets)) {
      if (verbose) {
        warning(
          "Only ", length(targets_intersect), " out of ", length(targets), " candidate regulators are in `matrix`."
        )
      }
    }
  }

  if (!is.null(regulators)) {
    regulators_intersect <- intersect(regulators, colnames(matrix))
    if (length(regulators_intersect) == 0) {
      stop("The input genes must contain at least 1 regulator.")
    }

    if (length(regulators_intersect) == 1) {
      if (verbose) {
        message("Only 1 regulator found in `matrix`, `weight` calculated by `cor`.")
      }
    }

    if (length(regulators_intersect) < length(regulators)) {
      if (verbose) {
        warning("Only ", length(regulators_intersect), " out of ", length(regulators), " candidate regulators are in the expression matrix.")
      }
    }
  }

  if (!is.numeric(cores) || cores < 1) {
    stop("`cores` should be a stricly positive integer.")
  }

  if (verbose) {
    message("Using `", penalty, "` penalty.")
    if (cross_validation) message("Using cross validation.")
  }
}

#' @title Switch weight table to matrix
#'
#' @param weight_table The weight data table of network.
#' @inheritParams net.format
#'
#' @return Weight matrix
#' @export
#'
#' @examples
#' library(inferCSN)
#' data("example_matrix")
#' weight_table <- inferCSN(example_matrix)
#' head(weight_table)
#'
#' table.to.matrix(weight_table)[1:6, 1:6]
#'
#' table.to.matrix(
#'   weight_table,
#'   regulators = c("g1", "g2"),
#'   targets = c("g3", "g4")
#' )
table.to.matrix <- function(
    weight_table,
    regulators = NULL,
    targets = NULL) {
  weight_table <- net.format(
    weight_table,
    abs_weight = FALSE
  )
  weight_matrix <- .Call(
    "_inferCSN_table_to_matrix",
    PACKAGE = "inferCSN",
    weight_table
  )
  weight_matrix <- filter_sort_matrix(
    weight_matrix,
    regulators = regulators,
    targets = targets
  )

  return(weight_matrix)
}

#' @title Filter and sort matrix
#'
#' @param weight_matrix The matrix of network weight.
#' @inheritParams net.format
#'
#' @return Filtered and sorted matrix
#' @export
#'
#' @examples
#' library(inferCSN)
#' data("example_matrix")
#' weight_table <- inferCSN(example_matrix)
#' weight_matrix <- table.to.matrix(weight_table)
#' filter_sort_matrix(weight_matrix)[1:6, 1:6]
#'
#' filter_sort_matrix(
#'   weight_matrix,
#'   regulators = c("g1", "g2"),
#'   targets = c("g3", "g4")
#' )
filter_sort_matrix <- function(
    weight_matrix,
    regulators = NULL,
    targets = NULL) {
  weight_matrix[is.na(weight_matrix)] <- 0
  if (!is.null(regulators)) {
    regulators <- intersect(rownames(weight_matrix), regulators)
  } else {
    regulators <- rownames(weight_matrix)
  }
  if (!is.null(targets)) {
    targets <- intersect(colnames(weight_matrix), targets)
  } else {
    targets <- colnames(weight_matrix)
  }

  unique_regulators <- gtools::mixedsort(unique(regulators))
  unique_targets <- gtools::mixedsort(unique(targets))
  weight_matrix <- weight_matrix[unique_regulators, unique_targets]

  return(weight_matrix)
}

#' @title Format weight table
#'
#' @param weight_table The weight data table of network.
#' @param regulators Regulators list.
#' @param targets Targets list.
#' @param abs_weight Logical value, whether to perform absolute value on weights,
#'  default set to `TRUE`, and when set `abs_weight` to `TRUE`,
#'  the output of weight table will create a new column named `Interaction`.
#'
#' @return Format weight table
#' @export
#'
#' @examples
#' library(inferCSN)
#' data("example_matrix")
#' weight_table <- inferCSN(example_matrix)
#'
#' net.format(
#'   weight_table,
#'   regulators = c("g1")
#' )
#'
#' net.format(
#'   weight_table,
#'   regulators = c("g1"),
#'   abs_weight = FALSE
#' )
#'
#' net.format(
#'   weight_table,
#'   targets = c("g3")
#' )
#'
#' net.format(
#'   weight_table,
#'   regulators = c("g1", "g3"),
#'   targets = c("g3", "g5")
#' )
net.format <- function(
    weight_table,
    regulators = NULL,
    targets = NULL,
    abs_weight = TRUE) {
  colnames(weight_table) <- c("regulator", "target", "weight")
  weight_table$weight <- as.numeric(weight_table$weight)
  weight_table <- dplyr::filter(weight_table, weight != 0)
  if (!is.null(regulators)) {
    weight_table <- purrr::map_dfr(
      unique(regulators), function(x) {
        dplyr::filter(weight_table, regulator == x)
      }
    )
  }
  if (!is.null(targets)) {
    weight_table <- purrr::map_dfr(
      unique(targets), function(x) {
        dplyr::filter(weight_table, target == x)
      }
    )
  }

  if (abs_weight) {
    weight_table$Interaction <- ifelse(
      weight_table$weight < 0, "Repression", "Activation"
    )
    weight_table$weight <- abs(weight_table$weight)
  }

  weight_table <- weight_table[order(
    abs(as.numeric(weight_table$weight)),
    decreasing = TRUE
  ), ]

  return(weight_table)
}

#' @title Extracts a specific solution in the regularization path
#'
#' @param object The output of model.fit or inferCSN.cvfit
#' @param lambda The value of lambda at which to extract the solution
#' @param gamma The value of gamma at which to extract the solution
#' @param supportSize The number of non-zeros each solution extracted will contain
#' @param ... Other parameters
#'
#' @method coef SRM_fit
#'
#' @return Return the specific solution
#' @export
coef.SRM_fit <- function(
    object,
    lambda = NULL,
    gamma = NULL,
    supportSize = NULL,
    ...) {
  if (!is.null(supportSize) && !is.null(lambda)) {
    stop("If 'supportSize' is provided to 'coef' only 'gamma' can also be provided.")
  }

  if (is.null(lambda) && is.null(gamma) && is.null(supportSize)) {
    # If all three are null, return all solutions
    t <- do.call(cbind, object$beta)
    if (object$settings$intercept) {
      intercepts <- unlist(object$a0)
      t <- rbind(intercepts, t)
    }
    return(t)
  }

  if (is.null(gamma)) gamma <- object$gamma[1]

  diffGamma <- abs(object$gamma - gamma)
  gammaindex <- which(diffGamma == min(diffGamma))

  indices <- NULL
  if (!is.null(lambda)) {
    diffLambda <- abs(lambda - object$lambda[[gammaindex]])
    indices <- which(diffLambda == min(diffLambda))
  } else if (!is.null(supportSize)) {
    diffSupportSize <- abs(supportSize - object$suppSize[[gammaindex]])
    indices <- which(diffSupportSize == min(diffSupportSize))
  } else {
    indices <- seq_along(object$lambda[[gammaindex]])
  }

  if (object$settings$intercept) {
    t <- rbind(
      object$a0[[gammaindex]][indices],
      object$beta[[gammaindex]][, indices, drop = FALSE]
    )
    rownames(t) <- c(
      "Intercept",
      paste0(rep("V", object$p), 1:object$p)
    )
  } else {
    t <- object$beta[[gammaindex]][, indices, drop = FALSE]
    rownames(t) <- paste0(rep("V", object$p), 1:object$p)
  }
  t
}

#' @rdname coef.SRM_fit
#'
#' @method coef SRM_fit_CV
#'
#' @return Return the specific solution
#' @export
coef.SRM_fit_CV <- function(
    object,
    lambda = NULL,
    gamma = NULL,
    ...) {
  coef.SRM_fit(object$fit, lambda, gamma, ...)
}

#' @title Prints a summary of model.fit
#'
#' @param x The output of model.fit or inferCSN.cvfit
#' @param ... Other parameters
#'
#' @method print SRM_fit
#'
#' @return Return information of model.fit
#' @export
print.SRM_fit <- function(x, ...) {
  gammas <- rep(x$gamma, times = lapply(x$lambda, length))
  data.frame(
    lambda = unlist(x["lambda"]),
    gamma = gammas,
    suppSize = unlist(x["suppSize"]),
    row.names = NULL
  )
}

#' @rdname print.SRM_fit
#'
#' @method print SRM_fit_CV
#'
#' @return Return information of model.fit
#' @export
print.SRM_fit_CV <- function(x, ...) {
  print.SRM_fit(x$fit)
}

#' @title Predict Response
#'
#' @description Predicts response for a given sample
#'
#' @param object The output of model.fit
#' @param newx A matrix on which predictions are made. The matrix should have p columns
#' @param lambda The value of lambda to use for prediction.
#' A summary of the lambdas in the regularization path can be obtained using \code{print(fit)}
#' @param gamma The value of gamma to use for prediction.
#' A summary of the gammas in the regularization path can be obtained using \code{print(fit)}
#' @param ... Other parameters
#'
#' @method predict SRM_fit
#'
#' @details
#' If both lambda and gamma are not supplied, then a matrix of predictions for all the solutions in the regularization path is returned.
#' If lambda is supplied but gamma is not, the smallest value of gamma is used.
#' In case of logistic regression, probability values are returned
#'
#' @return Return predict value
#' @export
predict.SRM_fit <- function(
    object,
    newx,
    lambda = NULL,
    gamma = NULL,
    ...) {
  beta <- coef.SRM_fit(object, lambda, gamma)
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
  prediction
}

#' @rdname predict.SRM_fit
#'
#' @method predict SRM_fit_CV
#'
#' @return Return the predict value
#' @export
predict.SRM_fit_CV <- function(
    object,
    newx,
    lambda = NULL,
    gamma = NULL,
    ...) {
  predict.SRM_fit(object$fit, newx, lambda, gamma, ...)
}

is.scalar <- function(x) {
  is.atomic(x) && length(x) == 1L && !is.character(x) && Im(x) == 0 && !is.nan(x) && !is.na(x)
}

#' normalization
#'
#' @param x A numeric vector.
#' @param method Method for normalization.
#'
#' @return Normalized vector
#' @export
normalization <- function(
    x,
    method = "max_min") {
  na_index <- which(is.na(x))
  x[na_index] <- 0
  x <- switch(
    EXPR = method,
    "max_min" = ((x - min(x)) / (max(x) - min(x))),
    "max" = (x / max(abs(x))),
    "sum" = (x / sum(abs(x))),
    "softmax" = .softmax(x)
  )
  x[na_index] <- NA

  return(x)
}

.softmax <- function(x) {
  abs_x <- abs(x)
  exp_abs_x <- exp(abs_x)
  sum_exp_abs_x <- sum(exp_abs_x)
  softmax_values <- exp_abs_x / sum_exp_abs_x
  result <- softmax_values * sign(x)

  return(result)
}
