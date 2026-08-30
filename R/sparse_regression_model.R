#' @title Fit greedy L0 regression
#' @description Fits a standardized sparse linear model by BIC-decreasing support updates.
#'
#' @param x Numeric predictor matrix with observations in rows.
#' @param y Numeric response vector.
#' @param max_support_size Optional positive integer support limit.
#' @param min_improvement Minimum RSS improvement considered.
#' @param verbose Logical verbosity flag.
#'
#' @return A list with `model`, `metrics`, and `coefficients`.
#' @export
fit_greedy_l0 <- function(
  x,
  y,
  max_support_size = NULL,
  min_improvement = 1e-10,
  verbose = TRUE
) {
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  y <- as.numeric(y)
  if (length(y) != nrow(x)) {
    stop("`y` must have one value for every row of `x`.", call. = FALSE)
  }
  if (is.null(colnames(x))) {
    colnames(x) <- paste0("V", seq_len(ncol(x)))
  }
  max_support_size <- validate_max_support_size(max_support_size)
  fit <- solve_greedy_l0(
    x = x,
    y = y,
    max_support_size = if (is.null(max_support_size)) 0L else max_support_size,
    min_improvement = as.numeric(min_improvement)
  )
  coefficient <- as.numeric(fit$coefficient)
  names(coefficient) <- colnames(x)
  support <- as.integer(fit$support)
  deletion_delta_bic <- rep(NA_real_, length(coefficient))
  if (length(support)) {
    standardized <- matrix(0, nrow = nrow(x), ncol = ncol(x))
    for (column in seq_len(ncol(x))) {
      values <- x[, column]
      finite <- is.finite(values)
      if (sum(finite) < 2L) next
      centered <- numeric(length(values))
      centered[finite] <- values[finite] - mean(values[finite])
      scale_value <- sqrt(sum(centered^2) / (sum(finite) - 1L))
      if (is.finite(scale_value) && scale_value > 0) {
        standardized[, column] <- centered / scale_value
      }
    }
    selected_gram <- crossprod(standardized[, support, drop = FALSE])
    inverse <- tryCatch(
      chol2inv(chol((selected_gram + t(selected_gram)) / 2)),
      error = function(error) NULL
    )
    if (is.null(inverse) || any(!is.finite(inverse)) ||
      any(diag(inverse) <= 0)) {
      stop("Selected support Gram matrix is not identifiable.", call. = FALSE)
    }
    beta <- coefficient[support]
    removed_rss <- as.numeric(fit$rss) + beta^2 / diag(inverse)
    removed_bic <- as.integer(fit$n_obs) * log(pmax(
      removed_rss / as.integer(fit$n_obs), 1e-12
    )) + (length(support) - 1L) * log(as.integer(fit$n_obs))
    delta <- removed_bic - as.numeric(fit$bic)
    tolerance <- 1e-8 * (1 + abs(as.numeric(fit$bic)))
    if (any(delta < -tolerance)) {
      stop("Selected support is not deletion-local-optimal.", call. = FALSE)
    }
    deletion_delta_bic[support] <- pmax(0, delta)
  }
  list(
    model = list(
      algorithm = "greedy_l0",
      support = colnames(x)[support],
      bic = as.numeric(fit$bic),
      rss = as.numeric(fit$rss),
      n_obs = as.integer(fit$n_obs)
    ),
    metrics = list(r_squared = as.numeric(fit$r_squared)),
    coefficients = list(
      variable = colnames(x),
      coefficient = coefficient,
      deletion_delta_bic = deletion_delta_bic
    )
  )
}

#' @title Fit batched greedy L0 regressions
#' @description Fits multiple responses from shared predictor cross-products.
#'
#' @param gram Numeric square predictor cross-product matrix.
#' @param xty Numeric predictor-by-response cross-product matrix.
#' @param response_ss Numeric response sums of squares.
#' @param candidates List of one-based predictor indices, one element per
#' response.
#' @param n_obs Number of observations used for the cross-products.
#' @param max_support_size Optional positive integer support limit.
#' @param min_improvement Minimum RSS reduction proposed during forward search.
#'
#' @return Batched coefficients, deletion evidence, and fit summaries.
#' @export
fit_greedy_l0_batch <- function(
  gram,
  xty,
  response_ss,
  candidates,
  n_obs,
  max_support_size = NULL,
  min_improvement = 1e-10
) {
  gram <- as.matrix(gram)
  xty <- as.matrix(xty)
  storage.mode(gram) <- storage.mode(xty) <- "double"
  max_support_size <- validate_max_support_size(max_support_size)
  solve_greedy_l0_batch(
    gram,
    xty,
    as.numeric(response_ss),
    candidates,
    as.integer(n_obs),
    if (is.null(max_support_size)) 0L else max_support_size,
    as.numeric(min_improvement)
  )
}

validate_max_support_size <- function(max_support_size) {
  if (is.null(max_support_size)) {
    return(NULL)
  }
  if (!is.numeric(max_support_size) || length(max_support_size) != 1L ||
    !is.finite(max_support_size) || max_support_size < 1 ||
    max_support_size != as.integer(max_support_size)) {
    stop(
      "`max_support_size` must be `NULL` or one positive integer.",
      call. = FALSE
    )
  }
  as.integer(max_support_size)
}
