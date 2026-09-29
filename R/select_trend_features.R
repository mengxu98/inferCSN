#' @title Select features with significant trends
#'
#' @md
#' @param x Numeric matrix with observations in rows and features in columns,
#'   or NULL when supplying `statistics`.
#' @param pseudotime Numeric coordinate per observation.
#' @param method Fitting method: `"pretsa"` or `"gam"`.
#' @param padjust_threshold Strict adjusted P-value cutoff between zero and one.
#' @param n_candidates Maximum number of selected features, or NULL for no cap.
#' @param statistics Precomputed table with feature row names and `padjust`
#'   (or `pvalue` if `padjust` is absent).
#' @param ... Numerical options passed to [thisutils::fit_trends()].
#'
#' @return A list with `features`, `statistics`, and `fit` (NULL for precomputed
#'   input). Nonfinite P-values are excluded. Input order is retained unless
#'   the cap is exceeded; then the smallest P-values are selected with stable ties.
#'
#' @export
#'
#' @examples
#' t <- seq(0, 1, length.out = 40)
#' x <- cbind(a = sin(t * 5), b = cos(t * 3))
#' select_trend_features(x, t)$features
select_trend_features <- function(
  x = NULL, pseudotime = NULL, method = c("pretsa", "gam"),
  padjust_threshold = 0.05, n_candidates = NULL, statistics = NULL, ...
) {
  method <- match.arg(method)
  if (length(padjust_threshold) != 1L || !is.finite(padjust_threshold) ||
    padjust_threshold < 0 || padjust_threshold > 1) {
    stop("`padjust_threshold` must be in [0, 1].", call. = FALSE)
  }
  if (!is.null(n_candidates) && (length(n_candidates) != 1L ||
    !is.finite(n_candidates) || n_candidates < 1 || n_candidates != floor(n_candidates))) {
    stop("`n_candidates` must be a positive integer or NULL.", call. = FALSE)
  }
  fit <- NULL
  if (is.null(statistics)) {
    fit <- thisutils::fit_trends(t(as.matrix(x)), pseudotime, method = method, ...)
    statistics <- fit$statistics
  } else if (length(list(...))) {
    stop("Numerical options cannot be used with precomputed statistics.", call. = FALSE)
  }
  if (!is.data.frame(statistics) || is.null(rownames(statistics)) || anyDuplicated(rownames(statistics))) {
    stop("`statistics` must have unique feature row names.", call. = FALSE)
  }
  p <- if ("padjust" %in% names(statistics)) statistics$padjust else statistics$pvalue
  if (!is.numeric(p) || length(p) != nrow(statistics)) {
    stop("`statistics` needs a numeric padjust or pvalue column.", call. = FALSE)
  }
  keep <- which(is.finite(p) & p < padjust_threshold)
  if (!is.null(n_candidates) && length(keep) > n_candidates) {
    keep <- keep[order(p[keep], seq_along(keep))][seq_len(n_candidates)]
  }
  list(features = rownames(statistics)[keep], statistics = statistics, fit = fit)
}
