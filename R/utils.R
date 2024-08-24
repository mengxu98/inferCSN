#' @title Print diagnostic message
#'
#' @param ... Text to print.
#' @param verbose Logical value, default is *`TRUE`*.
#' Whether to print the message.
#' @param message_type Type of message, default is *`info`*.
#' Could be choose one of *`info`*, *`warning`*, and *`error`*.
#' @param cli_model Logical value, default is *`TRUE`*.
#' Whether to use the `cli` package to print the message.
#' Add because the message is printed by \code{\link[base]{message}},
#' the message could be suppressed by \code{\link[base]{suppressMessages}}.
#'
#' @md
#' @export
#' @examples
#' log_message("Hello, ", "world!")
#' suppressMessages(log_message("Hello, ", "world!"))
#' log_message("Hello, world!", verbose = FALSE)
#' log_message("Hello, world!", verbose = TRUE, message_type = "warning")
log_message <- function(
    ...,
    verbose = TRUE,
    message_type = "info",
    cli_model = TRUE) {
  if (message_type == "error") {
    stop(...)
  }
  if (verbose) {
    if (cli_model) {
      switch(
        EXPR = message_type,
        "info" = cli::cli_alert_success(paste0(...)),
        "warning" = cli::cli_alert_warning(paste0("Warning: ", ...))
      )
    } else {
      switch(
        EXPR = message_type,
        "info" = message(paste0(...)),
        "warning" = message(paste0("Warning: ", ...))
      )
    }
  }
}

#' @title Parallelize a function
#'
#' @inheritParams inferCSN
#' @param x A vector or list to apply over.
#' @param fun The function to be applied to each element.
#' @param export_fun The functions to export the function to workers.
#'
#' @md
#' @return A list of computed results
#'
#' @export
parallelize_fun <- function(
    x,
    fun,
    cores = 1,
    export_fun = NULL,
    verbose = TRUE) {
  if (cores == 1) {
    log_message(
      "Using 1 core.",
      verbose = verbose
    )
    if (verbose) {
      return(pbapply::pblapply(X = x, FUN = fun))
    }
    if (!verbose) {
      return(base::lapply(X = x, FUN = fun))
    }
  }

  if (cores > 1) {
    doParallel::registerDoParallel(cores = cores)
    log_message(
      "Using ", foreach::getDoParWorkers(), " cores.",
      verbose = verbose
    )

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
  log_message(
    "Checking input parameters.",
    verbose = verbose
  )

  if (length(dim(matrix)) != 2) {
    log_message(
      "The input matrix must be a two-dimensional matrix.",
      message_type = "error",
      verbose = verbose
    )
  }

  if (is.null(colnames(matrix))) {
    log_message(
      "The input matrix must contain the names of the genes as colnames.",
      message_type = "error"
    )
  }

  match.arg(penalty, c("L0", "L0L1", "L0L2"))
  match.arg(algorithm, c("CD", "CDPSI"))

  if (!is.numeric(seed)) {
    seed <- 1
    log_message(
      "random seed is not a valid value, initialize it to 1.",
      message_type = "warning",
      verbose = verbose
    )
  }

  if (!(is.numeric(percent_samples) && percent_samples > 0 && percent_samples <= 1)) {
    log_message(
      "Please set `percent_samples` value between: (0, 1].",
      message_type = "error"
    )
  }

  if (!is.null(targets)) {
    intersect_targets <- intersect(targets, colnames(matrix))
    if (length(intersect_targets) == 0) {
      log_message(
        "The input genes must contain at least 1 target.",
        message_type = "error"
      )
    }

    if (length(intersect_targets) < length(targets)) {
      log_message(
        length(intersect_targets), " out of ",
        length(targets), " candidate targets are in the input matrix.",
        message_type = "warning",
        verbose = verbose
      )
    }
  }

  if (!is.null(regulators)) {
    intersect_regulators <- intersect(regulators, colnames(matrix))
    if (length(intersect_regulators) == 0) {
      log_message(
        "The input genes must contain at least 1 regulator.",
        message_type = "error"
      )
    }

    if (length(intersect_regulators) < length(regulators)) {
      log_message(
        length(intersect_regulators), " out of ",
        length(regulators), " candidate regulators are in the input matrix.",
        message_type = "warning",
        verbose = verbose
      )
    }
  }

  if (!is.numeric(cores) || cores < 1) {
    log_message(
      "`cores` should be a stricly positive integer.",
    )
  }

  log_message(
    "Using ", penalty, " sparse regression model.",
    verbose = verbose
  )
  if (cross_validation) {
    log_message(
      "Using cross validation, and setting ", n_folds, " folds.",
      verbose = verbose
    )
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
#' @param parallel Logical value, default is *`FALSE`*.
#' Setting to parallelize the computation with \code{\link[RcppParallel]{setThreadOptions}}.
#' @param sparse Logical value, default is *`FALSE`*, whether to output a sparse matrix.
#'
#' @md
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
#' \dontrun{
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
#' }
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

#' @title Format network table
#'
#' @param network_table The weight data table of network.
#' @param regulators Regulators list.
#' @param targets Targets list.
#' @param abs_weight Logical value, default is *`TRUE`*, whether to perform absolute value on weights,
#'  and when set `abs_weight` to *`TRUE`*,
#'  the output of weight table will create a new column named `Interaction`.
#'
#' @md
#' @return Formated network table
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

#' @title Filter and sort matrix
#'
#' @inheritParams network_format
#' @param network_matrix The matrix of network weight.
#'
#' @return Filtered and sorted matrix
#' @export
#'
#' @examples
#' data("example_matrix")
#' network_table <- inferCSN(example_matrix)
#' network_matrix <- table_to_matrix(network_table)
#' filter_sort_matrix(network_matrix)[1:6, 1:6]
#'
#' filter_sort_matrix(
#'   network_matrix,
#'   regulators = c("g1", "g2"),
#'   targets = c("g3", "g4")
#' )
filter_sort_matrix <- function(
  network_matrix,
  regulators = NULL,
  targets = NULL) {
network_matrix[is.na(network_matrix)] <- 0
if (!is.null(regulators)) {
  regulators <- intersect(rownames(network_matrix), regulators)
} else {
  regulators <- rownames(network_matrix)
}
if (!is.null(targets)) {
  targets <- intersect(colnames(network_matrix), targets)
} else {
  targets <- colnames(network_matrix)
}

unique_regulators <- gtools::mixedsort(unique(regulators))
unique_targets <- gtools::mixedsort(unique(targets))
network_matrix <- network_matrix[unique_regulators, unique_targets]

return(network_matrix)
}

#' @title Switch network table to matrix
#'
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
  network_matrix <- .Call(
    "_inferCSN_tableToMatrix",
    PACKAGE = "inferCSN",
    network_table
  )
  network_matrix <- filter_sort_matrix(
    network_matrix,
    regulators = regulators,
    targets = targets
  )

  return(network_matrix)
}

#' @title Extracts a specific solution in the regularization path
#'
#' @inheritParams inferCSN
#' @param object The output of \code{\link{fit_sparse_regression}}.
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

#' @title Prints a summary of `fit_sparse_regression`
#'
#' @param x The output of \code{\link{fit_sparse_regression}}.
#' @param ... Other parameters
#'
#' @method print srm
#'
#' @md
#' @return Return information of `fit_sparse_regression`
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
#' @param object The output of fit_sparse_regression.
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
#' If both lambda and gamma are not supplied, then a matrix of predictions for all the solutions in the regularization path is returned.
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
    "maximum" = {
      x / max(abs(x))
    },
    "sum" = {
      x / sum(abs(x))
    },
    "softmax" = {
      exp(x - max(x)) / sum(exp(x - max(x)))
    },
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

.is_scalar <- function(x) {
  is.atomic(x) && length(x) == 1L && !is.character(x) && Im(x) == 0 && !is.nan(x) && !is.na(x)
}

.rmse <- function(true, pred) {
  return(
    sqrt(mean((true - pred)^2))
  )
}

.sse <- function(y_true, y_pred) {
  return(
    sum((y_true - y_pred)**2)
  )
}

.rse <- function(y_true, y_pred) {
  return(
    .sse(y_true, y_pred) / .sse(y_true, mean(y_true))
  )
}

#' @title \eqn{R^2} (coefficient of determination)
#'
#' @param y_true A numeric vector with ground truth values.
#' @param y_pred A numeric vector with predicted values.
r_square <- function(y_true, y_pred) {
  return(
    1 - .rse(y_true, y_pred)
  )
}
