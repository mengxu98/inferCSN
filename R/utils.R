#' @title Apply function over a List or Vector
#'
#' @param x A vector or list to apply over.
#' @param fun The function to be applied to each element.
#' @param cores Number of CPU cores used.
#' Setting to parallelize the computation with \code{\link[foreach]{foreach}}.
#' @param export_fun export_fun.
#' @param verbose Logical. Whether to print progress bar.
#' Only works in sequential mode.
#'
#' @return A list of results.
#'
#' @export
parallelize_fun <- function(
    x,
    fun,
    cores = 1,
    export_fun = NULL,
    verbose = TRUE) {
  if (cores == 1 && verbose) {
    return(pbapply::pblapply(X = x, FUN = fun))
  }
  if (cores == 1 && !verbose) {
    return(base::lapply(X = x, FUN = fun))
  }
  if (cores > 1) {
    doParallel::registerDoParallel(cores = cores)
    if (verbose) {
      message("Using ", foreach::getDoParWorkers(), " cores.")
    }

    "%dopar%" <- foreach::"%dopar%"
    output_list <- foreach::foreach(
      i = 1:length(x),
      .export = export_fun
    ) %dopar% {
      fun(x[[i]])
    }
    names(output_list) <- names(x)

    doParallel::stopImplicitCluster()

    return(output_list)
  }
}

.check_parameters <- function(
    matrix,
    penalty,
    algorithm,
    cross_validation,
    seed,
    n_folds,
    percent_samples,
    r_threshold,
    regulators,
    targets,
    regulators_num,
    verbose,
    cores,
    ...) {
  if (verbose) {
    message("Checking input parameters.")
  }

  if (length(dim(matrix)) != 2) {
    stop(
      "The input matrix must be a two-dimensional matrix."
    )
  }

  if (is.null(colnames(matrix))) {
    stop(
      "The input matrix must contain the names of the genes as colnames."
    )
  }

  match.arg(penalty, c("L0", "L0L1", "L0L2"))
  match.arg(algorithm, c("CD", "CDPSI"))

  if (!is.numeric(seed)) {
    seed <- 1
    if (verbose) {
      message("Warning: random seed is not a valid value, initialize it to 1.")
    }
  }

  if (!(is.numeric(percent_samples) && percent_samples > 0 && percent_samples <= 1)) {
    stop("Please set `percent_samples` value between: (0, 1].")
  }

  if (!is.null(targets)) {
    intersect_targets <- intersect(targets, colnames(matrix))
    if (length(intersect_targets) == 0) {
      stop("The input genes must contain at least 1 target.")
    }

    if (length(intersect_targets) < length(targets)) {
      if (verbose) {
        message(
          "Warning: ",
          length(intersect_targets), " out of ",
          length(targets), " candidate targets are in the input matrix."
        )
      }
    }
  }

  if (!is.null(regulators)) {
    intersect_regulators <- intersect(regulators, colnames(matrix))
    if (length(intersect_regulators) == 0) {
      stop("The input genes must contain at least 1 regulator.")
    }

    if (length(intersect_regulators) < length(regulators)) {
      if (verbose) {
        message(
          "Warning: ",
          length(intersect_regulators), " out of ",
          length(regulators), " candidate regulators are in the input matrix."
        )
      }
    }
  }

  if (!is.numeric(cores) || cores < 1) {
    stop("`cores` should be a stricly positive integer.")
  }

  if (verbose) {
    message("Using `", penalty, "` penalty.")
    if (cross_validation) {
      message("Using cross validation.")
    }
  }
}

.cores_detect <- function(
    cores = 1,
    num_session = NULL) {
  if (is.null(num_session)) {
    return(1)
  } else {
    cores <- min(
      (parallel::detectCores(logical = FALSE) - 1), cores, num_session
    )

    return(cores)
  }
}

#' @title Convert dgCMatrix into a dense matrix
#'
#' @param x A matrix.
#' @param parallel Setting to parallelize the computation with \code{\link[RcppParallel]{setThreadOptions}}.
#' @param sparse Logical. Whether to output a sparse matrix.
#' @export
#'
#' @examples
#' dims_i <- 2000
#' dims_j <- 2000
#' sparse_matrix <- Matrix::sparseMatrix(
#'   i = sample(1:dims_i, 500),
#'   j = sample(1:dims_j, 500),
#'   x = rnorm(500),
#'   dims = c(dims_i, dims_j),
#'   dimnames = list(
#'     paste0("a", rep(1:dims_i)),
#'     paste0("b", rep(1:dims_j))
#'   )
#' )
#'
#' system.time(as.matrix(sparse_matrix))
#' system.time(as_matrix(sparse_matrix))
#' system.time(as_matrix(sparse_matrix, parallel = TRUE))
#'
#' identical(
#'   as.matrix(sparse_matrix),
#'   as_matrix(sparse_matrix)
#' )
#'
#' identical(
#'   as.matrix(sparse_matrix),
#'   as_matrix(sparse_matrix, parallel = TRUE)
#' )
#'
#' identical(
#'   sparse_matrix,
#'   as_matrix(as.matrix(sparse_matrix), sparse = TRUE)
#' )
#'
#' network_table_1 <- inferCSN(
#'   as_matrix(example_matrix, sparse = TRUE)
#' )
#' network_table_2 <- inferCSN(
#'   as(example_matrix, "sparseMatrix")
#' )
#'
#' plot_scatter(
#'   data.frame(
#'     network_table_1$weight,
#'     network_table_2$weight
#'   ),
#'   legend_position = "none"
#' )
as_matrix <- function(
    x,
    parallel = FALSE,
    sparse = FALSE) {
  if (!methods::is(x, "sparseMatrix")) {
    if (sparse) {
      non_zero <- which(x != 0, arr.ind = TRUE)
      return(
        Matrix::sparseMatrix(
          i = non_zero[, 1],
          j = non_zero[, 2],
          x = x[non_zero],
          dims = dim(x),
          dimnames = dimnames(x)
        )
      )
    } else {
      return(Matrix::as.matrix(x))
    }
  } else {
    row_pos <- x@i
    col_pos <- findInterval(seq_along(x@x) - 1, x@p[-1])
    if (parallel) {
      matrix <- .Call(
        "_inferCSN_asMatrixParallel",
        PACKAGE = "inferCSN",
        row_pos,
        col_pos,
        x@x,
        x@Dim[1],
        x@Dim[2]
      )
    } else {
      matrix <- .Call(
        "_inferCSN_asMatrix",
        PACKAGE = "inferCSN",
        row_pos,
        col_pos,
        x@x,
        x@Dim[1],
        x@Dim[2]
      )
    }

    attr(matrix, "dimnames") <- list(x@Dimnames[[1]], x@Dimnames[[2]])

    return(matrix)
  }
}

#' @title Check sparsity of matrix
#'
#' @param x A matrix.
#'
#' @return Sparsity of matrix
#' @export
check_sparsity <- function(x) {
  if (methods::is(x, "sparseMatrix")) {
    non_zero_count <- Matrix::nnzero(x)
    total_counts <- prod(dim(x))
  } else {
    non_zero_count <- sum(x != 0)
    total_counts <- length(x)
  }

  sparsity_ratio <- non_zero_count / total_counts

  sparsity <- 1 - sparsity_ratio

  return(sparsity)
}

#' @title Switch weight table to matrix
#'
#' @param network_table The weight data table of network.
#' @inheritParams network_format
#'
#' @return Weight matrix
#' @export
#'
#' @examples
#' data("example_matrix")
#' network_table <- inferCSN(example_matrix)
#' head(network_table)
#'
#' table_to_matrix(network_table)[1:6, 1:6]
#'
#' table_to_matrix(
#'   network_table,
#'   regulators = c("g1", "g2"),
#'   targets = c("g3", "g4")
#' )
table_to_matrix <- function(
    network_table,
    regulators = NULL,
    targets = NULL) {
  network_table <- network_format(
    network_table,
    abs_weight = FALSE
  )
  weight_matrix <- .Call(
    "_inferCSN_tableToMatrix",
    PACKAGE = "inferCSN",
    network_table
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
#' @inheritParams network_format
#'
#' @return Filtered and sorted matrix
#' @export
#'
#' @examples
#' data("example_matrix")
#' network_table <- inferCSN(example_matrix)
#' weight_matrix <- table_to_matrix(network_table)
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
#' @param network_table The weight data table of network.
#' @param regulators Regulators list.
#' @param targets Targets list.
#' @param abs_weight Logical value, whether to perform absolute value on weights,
#'  default set to `TRUE`, and when set `abs_weight` to `TRUE`,
#'  the output of weight table will create a new column named `Interaction`.
#'
#' @md
#' @return Format weight table
#' @export
#'
#' @examples
#' data("example_matrix")
#' network_table <- inferCSN(example_matrix)
#'
#' network_format(
#'   network_table,
#'   regulators = c("g1")
#' )
#'
#' network_format(
#'   network_table,
#'   regulators = c("g1"),
#'   abs_weight = FALSE
#' )
#'
#' network_format(
#'   network_table,
#'   targets = c("g3")
#' )
#'
#' network_format(
#'   network_table,
#'   regulators = c("g1", "g3"),
#'   targets = c("g3", "g5")
#' )
network_format <- function(
    network_table,
    regulators = NULL,
    targets = NULL,
    abs_weight = TRUE) {
  colnames(network_table) <- c("regulator", "target", "weight")
  network_table$weight <- as.numeric(network_table$weight)
  network_table <- dplyr::filter(network_table, weight != 0)
  if (!is.null(regulators)) {
    network_table <- purrr::map_dfr(
      unique(regulators), function(x) {
        dplyr::filter(network_table, regulator == x)
      }
    )
  }
  if (!is.null(targets)) {
    network_table <- purrr::map_dfr(
      unique(targets), function(x) {
        dplyr::filter(network_table, target == x)
      }
    )
  }

  if (abs_weight) {
    network_table$Interaction <- ifelse(
      network_table$weight < 0, "Repression", "Activation"
    )
    network_table$weight <- abs(network_table$weight)
  }

  network_table <- network_table[order(
    abs(as.numeric(network_table$weight)),
    decreasing = TRUE
  ), ]
  rownames(network_table) <- NULL

  return(network_table)
}

#' @title Extracts a specific solution in the regularization path
#'
#' @param object The output of \code{\link{fit_sparse_regression}}.
#' @param lambda The value of lambda at which to extract the solution.
#' @param gamma The value of gamma at which to extract the solution.
#' @param regulators_num The number of non-zeros each solution extracted will contain.
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
    regulators_num = NULL,
    ...) {
  if (!is.null(regulators_num) && !is.null(lambda)) {
    stop("If 'regulators_num' is provided to 'coef' only 'gamma' can also be provided.")
  }

  if (is.null(lambda) && is.null(gamma) && is.null(regulators_num)) {
    # If all three are null, return all solutions
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
    diffLambda <- abs(lambda - object$lambda[[gamma_index]])
    indices <- which(diffLambda == min(diffLambda))
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

#' @rdname coef.SRM_fit
#' @method coef SRM_fit_CV
#' @export
coef.SRM_fit_CV <- function(
    object,
    lambda = NULL,
    gamma = NULL,
    ...) {
  coef.SRM_fit(object$fit, lambda, gamma, ...)
}

#' @title Prints a summary of \code{fit_sparse_regression}
#'
#' @param x The output of \code{\link{fit_sparse_regression}}.
#' @param ... Other parameters
#'
#' @method print SRM_fit
#'
#' @return Return information of \code{fit_sparse_regression}
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
#' @method print SRM_fit_CV
#' @export
print.SRM_fit_CV <- function(x, ...) {
  print.SRM_fit(x$fit)
}

#' @title Predict Response
#'
#' @description Predicts response for a given sample
#'
#' @param object The output of fit_sparse_regression
#' @param newx A matrix on which predictions are made. The matrix should have p columns
#' @param lambda The value of lambda to use for prediction.
#' A summary of the lambdas in the regularization path can be obtained using \code{\link{print.SRM_fit}}.
#' @param gamma The value of gamma to use for prediction.
#' A summary of the gammas in the regularization path can be obtained using \code{\link{print.SRM_fit}}.
#' @param ... Other parameters
#'
#' @method predict SRM_fit
#'
#' @details
#' If both lambda and gamma are not supplied, then a matrix of predictions for all the solutions in the regularization path is returned.
#' If lambda is supplied but gamma is not, the smallest value of gamma is used.
#' In case of logistic regression, probability values are returned.
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

  return(prediction)
}

#' @rdname predict.SRM_fit
#' @method predict SRM_fit_CV
#' @export
predict.SRM_fit_CV <- function(
    object,
    newx,
    lambda = NULL,
    gamma = NULL,
    ...) {
  predict.SRM_fit(object$fit, newx, lambda, gamma, ...)
}

#' @title Normalize numeric vector
#'
#' @param x Input numeric vector.
#' @param method Method used for normalization.
#' @param na_rm Whether to remove `NA` values,
#' and if setting TRUE, using `0` instead.
#'
#' @md
#' @return Normalized numeric vector
#' @export
#'
#' @examples
#' nums <- c(runif(2), NA, -runif(2))
#' nums
#' normalization(nums, method = "max_min")
#' normalization(nums, method = "maximum")
#' normalization(nums, method = "sum")
#' normalization(nums, method = "softmax")
#' normalization(nums, method = "z_score")
#' normalization(nums, method = "mad")
#' normalization(nums, method = "unit_vector")
#' normalization(nums, method = "unit_vector", na_rm = FALSE)
normalization <- function(
    x,
    method = "max_min",
    na_rm = TRUE) {
  method <- match.arg(
    method,
    c("max_min", "maximum", "sum", "softmax", "z_score", "mad", "unit_vector")
  )
  na_index <- which(is.na(x))
  x[na_index] <- 0
  x <- switch(
    EXPR = method,
    "max_min" = {
      (x - min(x)) / (max(x) - min(x))
    },
    "maximum" = (x / max(abs(x))),
    "sum" = (x / sum(abs(x))),
    "softmax" = .softmax(x),
    "z_score" = {
      (x - mean(x)) / stats::sd(x)
    },
    "mad" = {
      (x - stats::median(x)) / stats::mad(x)
    },
    "unit_vector" = {
      x / sqrt(sum(x^2))
    }
  )

  if (!na_rm) {
    x[na_index] <- NA
  }

  return(x)
}

.softmax <- function(x) {
  abs_x <- abs(x)
  softmax_values <- exp(abs_x) / sum(exp(abs_x))
  result <- softmax_values * sign(x)

  return(result)
}

.is_scalar <- function(x) {
  is.atomic(x) && length(x) == 1L && !is.character(x) && Im(x) == 0 && !is.nan(x) && !is.na(x)
}

.rmse <- function(true, pred) {
  sqrt(mean((true - pred)^2))
}

#' @title Sum of Squared Errors
#'
#' @param y_true A numeric vector with ground truth values.
#' @param y_pred A numeric vector with predicted values.
sse <- function(y_true, y_pred) {
  return(sum((y_true - y_pred)**2))
}

#' @title Relative Squared Error
#'
#' @inheritParams sse
rse <- function(y_true, y_pred) {
  return(sse(y_true, y_pred) / sse(y_true, mean(y_true)))
}

#' @title \eqn{R^2} (coefficient of determination)
#'
#' @inheritParams sse
r_square <- function(y_true, y_pred) {
  1 - rse(y_true, y_pred)
}
