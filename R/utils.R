utils::globalVariables(c(
  "x",
  "y",
  "xend",
  "yend",
  "regulator",
  "target",
  "weight",
  "Interaction",
  "name",
  "degree",
  "edges",
  "curvetype"
))

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
  if (verbose) message("Checking input parameters.")

  if (!is.matrix(matrix) && !is.array(matrix) && (length(dim(matrix)) != 2)) {
    stop(
      "Parameter matrix must be a two-dimensional matrix,
      where each column corresponds to a gene and each row corresponds to a sample/cell."
    )
  }

  if (is.null(colnames(matrix))) {
    stop("Parameter matrix must contain the names of the genes as colnames.")
  }

  # Check the penalty term of the regression model
  if (!any(c("L0", "L0L2") == penalty)) {
    stop(
      "inferCSN does not support ", penalty, " penalty regression.\n",
      "Please set penalty item as 'L0' or 'L0L2'."
    )
  }

  # Check the algorithm of the regression model
  if (!any(c("CD", "CDPSI") == algorithm)) {
    stop(
      "inferCSN does not support ", algorithm, " algorithmn.\n",
      "Please set algorithm item as 'CD' or 'CDPSI'."
    )
  }

  if (!is.numeric(seed)) {
    seed <- 1
    if (verbose) warning("Supplied seed is not a valid integer, initialize 'seed' to 1.")
  }

  if (!is.null(k_folds)) {
    if (!(k_folds > 0 && k_folds < 10)) {
      stop("Please set 'k_folds' value between: (0, 10).")
    }
  }

  if (!is.null(targets)) {
    if (!is.vector(targets)) {
      stop("Parameter 'targets' must a vector (of indices or gene names).")
    }

    if (is.numeric(targets)) {
      if (max(targets) > nrow(matrix)) stop("At least one index in 'targets' exceeds the number of genes.")
      if (min(targets) <= 0) stop("The indexes in 'targets' should be >=1.")
    }

    if (any(table(targets) > 1)) stop("Please provide each target only once.")

    if (is.character(targets)) {
      targetsInMatrix <- intersect(targets, colnames(matrix))
      if (length(targetsInMatrix) == 0) {
        stop("The genes must contain at least one target.")
      }

      if (length(targetsInMatrix) < length(targets)) {
        if (verbose) warning("Only ", length(targetsInMatrix), " out of ", length(targets), " candidate regulators are in the expression matrix.")
      }
    }
  }

  if (!is.null(regulators)) {
    if (is.list(regulators)) {
      regulators <- unique(unlist(regulators))
    }
    if (!is.null(regulators)) {
      if (length(regulators) < 2) stop("Provide at least 2 potential regulators.")

      if (!is.vector(regulators)) {
        stop("Parameter 'regulators' must a vector of indices or gene names.")
      }

      if (is.numeric(regulators)) {
        if (max(regulators) > nrow(matrix)) {
          stop("At least one index in 'regulators' exceeds the number of genes.")
        }
        if (min(regulators) <= 0) stop("The indexes in 'regulators' should be >=1.")
      }

      if (any(table(regulators) > 1)) stop("Please provide each regulator only once.")

      if (is.character(regulators)) {
        regulators_in_matrix <- intersect(regulators, colnames(matrix))
        if (length(regulators_in_matrix) < 2) {
          stop("Fewer than 2 regulators in the columns of expression matrix.")
        }

        if (length(regulators_in_matrix) < length(regulators)) {
          if (verbose) {
            warning("Only ", length(regulators_in_matrix), " out of ", length(regulators), " candidate regulators are in the expression matrix.")
          }
        }
      }
    }
  }

  if (!is.numeric(cores) || cores < 1) {
    stop(x = "Parameter cores should be a stricly positive integer.")
  }

  if (verbose) {
    message("All parameters check done.")
    message("Using '", penalty, "' penalty.")
    if (cross_validation) message("Using cross validation.")
  }
}

#' @title Switch weight table to matrix
#'
#' @param weight_table The weight data table of network.
#'
#' @return Weight matrix
#' @export
table.to.matrix <- function(weight_table) {
  .Call(
    "_inferCSN_table_to_matrix",
    PACKAGE = "inferCSN",
    weight_table
  )
}

#' @title Format weight table
#'
#' @param weight_table The weight data table of network.
#' @param regulators Regulators list.
#' @param targets Targets list.
#' @param abs_weight Logical value, whether to perform absolute value on weights.
#'
#' @return Format weight table
#' @export
#'
#' @examples
#' library(inferCSN)
#' data("example_matrix")
#' data("example_ground_truth")
#' weight_table <- inferCSN(example_matrix)
#'
#' net.format(
#'   weight_table,
#'   regulators = c("g1")
#' )
#' net.format(
#'   weight_table,
#'   targets = c("g3")
#' )
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
      })
  }
  if (!is.null(targets)) {
    weight_table <- purrr::map_dfr(
      unique(targets), function(x) {
        dplyr::filter(weight_table, target == x)
      })
  }
  weight_table$Interaction <- ifelse(weight_table$weight < 0, "Repression", "Activation")

  if (abs_weight) {
    weight_table$weight <- abs(weight_table$weight)
  }

  return(weight_table)
}

#' @title Extracts a specific solution in the regularization path
#'
#' @param object The output of inferCSN.fit or inferCSN.cvfit
#' @param lambda The value of lambda at which to extract the solution
#' @param gamma The value of gamma at which to extract the solution
#' @param supportSize The number of non-zeros each solution extracted will contain
#' @param ... Other parameters
#'
#' @method coef inferCSN
#'
#' @return Return the specific solution
#' @export
coef.inferCSN <- function(
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

#' @rdname coef.inferCSN
#'
#' @method coef inferCSNCV
#'
#' @return Return the specific solution
#' @export
coef.inferCSNCV <- function(
    object,
    lambda = NULL,
    gamma = NULL,
    ...) {
  coef.inferCSN(object$fit, lambda, gamma, ...)
}

#' @title Prints a summary of inferCSN.fit
#'
#' @param x The output of inferCSN.fit or inferCSN.cvfit
#' @param ... Other parameters
#'
#' @method print inferCSN
#'
#' @return Return the information of inferCSN.fit
#' @export
print.inferCSN <- function(x, ...) {
  gammas <- rep(x$gamma, times = lapply(x$lambda, length))
  data.frame(
    lambda = unlist(x["lambda"]),
    gamma = gammas,
    suppSize = unlist(x["suppSize"]),
    row.names = NULL
  )
}

#' @rdname print.inferCSN
#'
#' @method print inferCSNCV
#'
#' @return Return the information of inferCSN.fit
#' @export
print.inferCSNCV <- function(x, ...) {
  print.inferCSN(x$fit)
}

#' @title Predict Response
#'
#' @description Predicts the response for a given sample
#'
#' @param object The output of inferCSN.fit
#' @param newx A matrix on which predictions are made. The matrix should have p columns
#' @param lambda The value of lambda to use for prediction.
#' A summary of the lambdas in the regularization path can be obtained using \code{print(fit)}
#' @param gamma The value of gamma to use for prediction.
#' A summary of the gammas in the regularization path can be obtained using \code{print(fit)}
#' @param ... Other parameters
#'
#' @method predict inferCSN
#'
#' @details
#' If both lambda and gamma are not supplied, then a matrix of predictions for all the solutions in the regularization path is returned.
#' If lambda is supplied but gamma is not, the smallest value of gamma is used.
#' In case of logistic regression, probability values are returned
#'
#' @return Return the predict value
#' @export
predict.inferCSN <- function(
    object,
    newx,
    lambda = NULL,
    gamma = NULL,
    ...) {
  beta <- coef.inferCSN(object, lambda, gamma)
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

#' @rdname predict.inferCSN
#'
#' @method predict inferCSNCV
#'
#' @return Return the predict value
#' @export
predict.inferCSNCV <- function(
    object,
    newx,
    lambda = NULL,
    gamma = NULL,
    ...) {
  predict.inferCSN(object$fit, newx, lambda, gamma, ...)
}

is.scalar <- function(x) {
  is.atomic(x) && length(x) == 1L && !is.character(x) && Im(x) == 0 && !is.nan(x) && !is.na(x)
}
