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
    thisutils::log_message("`y` must have one value for every row of `x`.", message_type = "error")
  }
  if (is.null(colnames(x))) {
    colnames(x) <- paste0("V", seq_len(ncol(x)))
  }
  max_support_size <- validate_max_support_size(max_support_size)
  fit <- solve_greedy_l0(
    x = x,
    y = y,
    max_support_size = max_support_size,
    min_improvement = as.numeric(min_improvement)
  )
  coefficient <- as.numeric(fit$coefficient)
  names(coefficient) <- colnames(x)
  support <- as.integer(fit$support)
  deletion_delta_bic <- as.numeric(fit$deletion_delta_bic)
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
    max_support_size,
    as.numeric(min_improvement)
  )
}

validate_max_support_size <- function(max_support_size) {
  if (is.null(max_support_size)) {
    return(0L)
  }
  if (!is.numeric(max_support_size) || length(max_support_size) != 1L ||
    !is.finite(max_support_size) || max_support_size < 1 ||
    max_support_size != as.integer(max_support_size)) {
    thisutils::log_message(
      "`max_support_size` must be `NULL` or one positive integer.",
      message_type = "error"
    )
  }
  as.integer(max_support_size)
}
