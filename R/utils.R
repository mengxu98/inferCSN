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

#' @title Value selection operator
#'
#' @description
#' This operator returns the left side if it's not NULL,
#' otherwise it returns the right side.
#'
#' @param a The left side value to check
#' @param b The right side value to use if a is NULL
#'
#' @export
#'
#' @examples
#' NULL %s% 10
#' 5 %s% 10
`%s%` <- function(a, b) {
  if (is.null(a)) {
    return(b)
  } else {
    return(a)
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
      i = seq_along(x),
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
    subsampling_method,
    subsampling_ratio,
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
      message_type = "error"
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
      "initialize random seed to 1.",
      message_type = "warning",
      verbose = verbose
    )
  }

  match.arg(
    subsampling_method,
    c("sample", "meta_cells")
  )

  if (!(is.numeric(subsampling_ratio) && subsampling_ratio > 0 && subsampling_ratio <= 1)) {
    log_message(
      "Please set 'subsampling_ratio' value between: (0, 1].",
      message_type = "error"
    )
  }

  if (!is.null(regulators)) {
    intersect_regulators <- intersect(regulators, colnames(matrix))
    if (length(intersect_regulators) < 2) {
      log_message(
        "The input genes must contain at least 2 regulator.",
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
    } else {
      log_message(
        "Using ", length(intersect_regulators), " regulator(s).",
        verbose = verbose
      )
    }
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
    } else {
      log_message(
        "Using ", length(intersect_targets), " target(s).",
        verbose = verbose
      )
    }
  }

  if (!is.numeric(cores) || cores < 1) {
    log_message(
      "`cores` should be a positive integer, initialize it to 1.",
      message_type = "warning",
      verbose = verbose
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

#' @title Generate a simulated sparse matrix for single-cell data testing
#'
#' @param nrow Number of rows (genes) in the matrix.
#' @param ncol Number of columns (cells) in the matrix.
#' @param density Density of non-zero elements (default: 0.1, representing 90 sparsity).
#' @param distribution_fun Function to generate non-zero values.
#' @param seed Random seed for reproducibility.
#'
#' @return A sparse matrix of class "dgCMatrix"
#' @export
#'
#' @examples
#' simulate_sparse_matrix(2000, 500) |>
#'   check_sparsity()
simulate_sparse_matrix <- function(
    nrow,
    ncol,
    density = 0.1,
    distribution_fun = function(n) stats::rpois(n, lambda = 0.5) + 1,
    seed = 1) {
  set.seed(seed)

  nnz <- round(nrow * ncol * density)

  i <- sample(1:nrow, nnz, replace = TRUE)
  j <- sample(1:ncol, nnz, replace = TRUE)
  x <- distribution_fun(nnz)

  Matrix::sparseMatrix(
    i = i,
    j = j,
    x = x,
    dims = c(nrow, ncol),
    dimnames = list(
      paste0("cell_", 1:nrow),
      paste0("gene_", 1:ncol)
    )
  )
}

#' @title Correlation and covariance calculation for sparse matrix
#'
#' @inheritParams sparse_cor
pearson_correlation <- function(x, y = NULL) {
  if (!methods::is(x, "sparseMatrix")) {
    stop("x should be a sparse matrix.")
  }
  if (!is.null(y) && !methods::is(y, "sparseMatrix")) {
    stop("y should be a sparse matrix.")
  }

  result <- sparseCovCor(x, y)

  return(
    list(
      cov = result$cov,
      cor = result$cor
    )
  )
}

#' @title Safe correlation function which returns a sparse matrix without missing values
#'
#' @param x Sparse matrix or character vector.
#' @param y Sparse matrix or character vector.
#' @param method Method to use for calculating the correlation coefficient.
#' @param allow_neg Logical. Whether to allow negative values or set them to 0.
#' @param remove_na Logical. Whether to replace NA values with 0.
#' @param remove_inf Logical. Whether to replace infinite values with 1.
#' @param ... Other arguments passed to \code{\link[stats]{cor}} function.
#'
#' @return A correlation matrix.
#'
#' @export
#'
#' @examples
#' m1 <- simulate_sparse_matrix(
#'   1000, 1000,
#'   density = 0.01
#' )
#' m2 <- simulate_sparse_matrix(
#'   1000, 500,
#'   density = 0.01
#' )
#'
#' all.equal(
#'   as.matrix(sparse_cor(m1)),
#'   cor(as_matrix(m1))
#' )
#' all.equal(
#'   as.matrix(sparse_cor(m1, m2)),
#'   cor(as_matrix(m1), as_matrix(m2))
#' )
#'
#' system.time(
#'   sparse_cor(m1)
#' )
#' system.time(
#'   cor(as_matrix(m1))
#' )
#' system.time(
#'   sparse_cor(m1, m2)
#' )
#' system.time(
#'   cor(as_matrix(m1), as_matrix(m2))
#' )
#'
#' # add missing values
#' m1[sample(1:500, 10)] <- NA
#' m2[sample(1:500, 10)] <- NA
#'
#' sparse_cor(m1, m2)[1:5, 1:5]
sparse_cor <- function(
    x,
    y = NULL,
    method = "pearson",
    allow_neg = TRUE,
    remove_na = TRUE,
    remove_inf = TRUE,
    ...) {
  if (!methods::is(x, "sparseMatrix")) {
    x <- as_matrix(x, sparse = TRUE)
  }

  if (!is.null(y)) {
    if (!methods::is(y, "sparseMatrix")) {
      y <- as_matrix(y, sparse = TRUE)
    }
    if (nrow(x) != nrow(y)) {
      stop("x and y must have the same number of rows.")
    }
  }

  corr_mat <- switch(
    EXPR = method,
    "pearson" = pearson_correlation(x, y)$cor,
    "spearman" = {
      if (is.null(y)) {
        stats::cor(
          as_matrix(x),
          method = "spearman",
          ...
        )
      } else {
        stats::cor(
          as_matrix(x),
          as_matrix(y),
          method = "spearman",
          ...
        )
      }
    },
    "kendall" = {
      if (is.null(y)) {
        stats::cor(
          as_matrix(x),
          method = "kendall",
          ...
        )
      } else {
        stats::cor(
          as_matrix(x),
          as_matrix(y),
          method = "kendall",
          ...
        )
      }
    }
  )

  if (is.null(y)) {
    colnames(corr_mat) <- colnames(x)
    rownames(corr_mat) <- colnames(x)
  } else {
    colnames(corr_mat) <- colnames(y)
    rownames(corr_mat) <- colnames(x)
  }

  if (remove_na) {
    corr_mat[is.na(corr_mat)] <- 0
  }
  if (remove_inf) {
    corr_mat[is.infinite(corr_mat)] <- 1
  }

  corr_mat <- as_matrix(corr_mat, sparse = TRUE)

  if (!allow_neg) {
    corr_mat[corr_mat < 0] <- 0
  }

  return(corr_mat)
}

#' @title Convert sparse matrix into dense matrix
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
#' sparse_matrix <- simulate_sparse_matrix(
#'   2000,
#'   2000,
#'   density = 0.01
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
#' network_table_0 <- inferCSN(example_matrix)
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
#'     network_table_0$weight,
#'     network_table_1$weight
#'   ),
#'   legend_position = "none"
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
      return(
        Matrix::Matrix(
          x,
          sparse = TRUE,
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
      matrix <- asMatrixParallel(
        row_pos,
        col_pos,
        x@x,
        x@Dim[1],
        x@Dim[2]
      )
    } else {
      matrix <- asMatrix(
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

#' @title Extracts a specific solution in the regularization path
#'
#' @inheritParams single_network
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

#' @title Normalize numeric vector
#'
#' @param x Input numeric vector.
#' @param method Method used for normalization.
#' @param na_rm Whether to remove `NA` values,
#' and if setting TRUE, using `0` instead.
#' @param ... Parameters for other methods.
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
    na_rm = TRUE,
    ...) {
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
