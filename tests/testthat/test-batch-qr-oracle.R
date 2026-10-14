test_that("batch supports and deletion evidence agree with direct QR refits", {
  for (seed in 3101:3108) {
    set.seed(seed)
    n <- 70L
    x <- matrix(rnorm(n * 7), nrow = n)
    x[, 6] <- x[, 1] + x[, 2] + rnorm(n, sd = 0.15)
    x[, 7] <- x[, 3] + rnorm(n, sd = 0.02)
    x <- scale(x)
    y <- scale(cbind(
      1.2 * x[, 1] - 0.7 * x[, 2] + rnorm(n, sd = 0.3),
      x[, 3] + 0.5 * x[, 4] + rnorm(n, sd = 0.2)
    ))
    candidates <- list(seq_len(7), c(1L, 3L, 4L, 5L, 7L))
    cap <- 4L
    fit <- fit_greedy_l0_batch(crossprod(x), crossprod(x, y), colSums(y^2),
      candidates,
      n_obs = n, max_support_size = cap
    )
    for (target in 1:2) {
      keep <- which(fit$target_index == target)
      selected <- fit$predictor_index[keep]
      qr_fit <- function(s) {
        if (!length(s)) {
          return(list(rss = sum(y[, target]^2), coefficient = numeric()))
        }
        model <- lm.fit(x[, s, drop = FALSE], y[, target])
        if (model$rank != length(s)) {
          return(list(rss = Inf, coefficient = numeric()))
        }
        list(rss = sum(model$residuals^2), coefficient = model$coefficients)
      }
      bic <- function(s) n * log(max(qr_fit(s)$rss / n, 1e-12)) + length(s) * log(n)
      expect_true(all(selected %in% candidates[[target]]))
      expect_lte(length(selected), cap)
      expect_equal(fit$rss[target], qr_fit(selected)$rss, tolerance = 1e-7)
      expect_equal(fit$bic[target], bic(selected), tolerance = 1e-7)
      expect_equal(fit$standardized_beta[keep], unname(qr_fit(selected)$coefficient), tolerance = 1e-7)
      for (j in seq_along(selected)) {
        delta <- bic(selected[-j]) - bic(selected)
        expect_gte(delta, -1e-7)
        expect_equal(fit$deletion_delta_bic[keep[j]], delta, tolerance = 1e-6)
      }
      outside <- setdiff(candidates[[target]], selected)
      neighbours <- lapply(seq_along(selected), function(j) selected[-j])
      if (length(selected) < cap) neighbours <- c(neighbours, lapply(outside, function(j) c(selected, j)))
      for (j in seq_along(selected)) neighbours <- c(neighbours, lapply(outside, function(k) c(selected[-j], k)))
      expect_true(all(vapply(neighbours, bic, numeric(1)) >= bic(selected) - 1e-6))
    }
  }
})

test_that("pivoted Gram fits restore predictor columns and deletion evidence", {
  gram <- matrix(c(100, 200, 0, 200, 900, 100, 0, 100, 400), 3)
  beta <- c(3, -1, 2)
  xty <- gram %*% beta
  y_ss <- sum(beta * xty) + 1
  expected_evidence <- 100 * log(1 + beta^2 / diag(solve(gram))) - log(100)

  for (candidate_order in list(1:3, 3:1, c(2L, 1L, 3L))) {
    fit <- fit_greedy_l0_batch(gram, xty, y_ss, list(candidate_order), n_obs = 100)
    expect_setequal(fit$predictor_index, 1:3)
    index <- match(1:3, fit$predictor_index)
    expect_equal(fit$standardized_beta[index], beta, tolerance = 1e-12)
    expect_equal(fit$rss, 1, tolerance = 1e-10)
    expect_equal(fit$deletion_delta_bic[index], expected_evidence, tolerance = 1e-10)
  }
})

test_that("constant, singleton, duplicate and wide predictors obey support contracts", {
  set.seed(3190)
  signal <- rnorm(12)
  fixtures <- list(
    matrix(1, 12, 2), matrix(signal, ncol = 1),
    cbind(signal, signal, rnorm(12)), matrix(rnorm(12 * 20), 12, 20)
  )
  for (x in fixtures) {
    fit <- fit_greedy_l0(x, signal + rnorm(12, sd = 0.1), verbose = FALSE)
    expect_lte(length(fit$model$support), min(ncol(x), nrow(x) - 2L))
    expect_true(all(is.finite(fit$coefficients$coefficient)))
    selected <- fit$coefficients$coefficient != 0
    expect_true(all(fit$coefficients$deletion_delta_bic[selected] >= 0))
  }
  empty <- inferCSN(cbind(a = rep(1, 20), b = rep(2, 20)), verbose = FALSE)
  expect_equal(nrow(empty), 0L)
})

test_that("long forward paths preserve coefficients and deletion evidence with pivoting", {
  for (seed in 3193:3196) {
    set.seed(seed)
    x <- matrix(rnorm(40 * 12), 40, 12)
    x[, 12] <- x[, 1] + x[, 2] + rnorm(40, sd = 0.1)
    x <- sweep(x, 2, rep(c(1, 4, 2), 4), "*")
    gram <- crossprod(x)
    beta <- rep(c(3, -2, 1), 4)
    xty <- gram %*% beta
    y_ss <- sum(beta * xty) + 1
    evidence <- 1000 * log(1 + beta^2 / diag(solve(gram))) - log(1000)
    for (candidate_order in list(1:12, 12:1)) {
      fit <- fit_greedy_l0_batch(gram, xty, y_ss, list(candidate_order), n_obs = 1000)
      expect_setequal(fit$predictor_index, 1:12)
      index <- match(1:12, fit$predictor_index)
      expect_equal(fit$standardized_beta[index], beta, tolerance = 1e-8)
      expect_equal(fit$rss, 1, tolerance = 1e-7)
      expect_equal(fit$deletion_delta_bic[index], evidence, tolerance = 1e-7)
    }
  }
})

test_that("large selected supports preserve coefficients and deletion evidence", {
  for (p in c(63L, 64L, 65L, 72L)) {
    set.seed(3197 + p)
    x <- matrix(rnorm(100 * p), 100, p)
    x <- sweep(x, 2, rep_len(c(1, 4, 2), p), "*")
    gram <- crossprod(x)
    beta <- rep_len(c(3, -2, 1), p)
    xty <- gram %*% beta
    y_ss <- sum(beta * xty) + 1
    evidence <- 1000 * log(1 + beta^2 / diag(solve(gram))) - log(1000)
    fit <- fit_greedy_l0_batch(
      gram, xty, y_ss, list(seq_len(p)),
      n_obs = 1000
    )
    expect_setequal(fit$predictor_index, seq_len(p))
    index <- match(seq_len(p), fit$predictor_index)
    expect_equal(fit$standardized_beta[index], beta, tolerance = 1e-8)
    expect_equal(fit$rss, 1, tolerance = 1e-7)
    expect_equal(fit$deletion_delta_bic[index], evidence, tolerance = 1e-7)
  }
})

test_that("static network output is invariant to core count and row ordering", {
  set.seed(3191)
  x <- matrix(rnorm(80 * 5), 80, 5, dimnames = list(paste0("c", 1:80), paste0("g", 1:5)))
  x[, 4] <- x[, 1] - x[, 2] + rnorm(80, sd = 0.1)
  one <- expect_network_across_cores(x)
  canonical <- function(z) {
    z <- z[order(z$regulator, z$target), , drop = FALSE]
    rownames(z) <- NULL
    z
  }
  expect_equal(canonical(one), canonical(inferCSN(x[80:1, ], cores = 1, verbose = FALSE)), tolerance = 1e-12)
})

test_that("batch agrees with a QR backward reference removing multiple redundant predictors", {
  set.seed(3192)
  n <- 100L
  x <- scale(matrix(rnorm(n * 6), n, 6))
  x[, 6] <- scale(x[, 1] + x[, 2] + rnorm(n, sd = .4))
  y <- as.numeric(scale(1.5 * x[, 1] - x[, 2] + rnorm(n, sd = .01)))
  bic <- function(s) {
    rss <- if (length(s)) sum(lm.fit(x[, s, drop = FALSE], y)$residuals^2) else sum(y^2)
    n * log(max(rss / n, 1e-12)) + length(s) * log(n)
  }
  support <- 1:6
  trajectory <- bic(support)
  repeat {
    trials <- lapply(seq_along(support), function(i) support[-i])
    scores <- vapply(trials, bic, numeric(1))
    if (!length(scores) || min(scores) >= tail(trajectory, 1) - 1e-8) break
    support <- trials[[which.min(scores)]]
    trajectory <- c(trajectory, min(scores))
  }
  expect_gte(length(trajectory) - 1L, 2L)
  expect_true(all(diff(trajectory) < 0))
  fit <- fit_greedy_l0_batch(crossprod(x), matrix(crossprod(x, y), ncol = 1), sum(y^2), list(1:6), n)
  expect_setequal(fit$predictor_index, support)
  expect_equal(fit$bic, tail(trajectory, 1), tolerance = 1e-7)
})

test_that("batch workspace reuse preserves changing candidate sets and target order", {
  set.seed(1013)
  x <- scale(matrix(rnorm(80 * 70), 80))
  y <- cbind(
    x[, 2] - x[, 8] + rnorm(80, sd = 0.2), rnorm(80), 0,
    x[, 5] + rnorm(80, sd = 0.3), rnorm(80)
  )
  gram <- crossprod(x)
  xty <- crossprod(x, y)
  ss <- colSums(y^2)
  candidates <- list(70:1, c(68L, 2L), 1:70, integer(), c(11L, 5L, 3L))
  batch <- fit_greedy_l0_batch(gram, xty, ss, candidates, nrow(x))
  for (i in seq_along(candidates)) {
    single <- fit_greedy_l0_batch(
      gram, xty[, i, drop = FALSE], ss[i],
      candidates[i], nrow(x)
    )
    edges <- batch$target_index == i
    for (field in c("predictor_index", "standardized_beta", "deletion_delta_bic")) {
      expect_identical(batch[[field]][edges], single[[field]])
    }
    for (field in c("candidate_size", "support_size", "rss", "bic", "status")) {
      expect_identical(batch[[field]][i], single[[field]])
    }
  }
})

test_that("batch validation retains selected-domain finite and index checks", {
  gram <- diag(3)
  xty <- matrix(c(1, 2, 3), ncol = 1)
  fit <- function(g = gram, c = list(1:3), cross = xty) {
    fit_greedy_l0_batch(g, cross, 10, c, 20)
  }
  for (bad in list(c(1L, 1L), c(1L, NA_integer_), 0L, -1L, 4L)) {
    expect_error(fit(c = list(bad)), "unique and one-based")
  }
  for (bad in c(NA_real_, NaN, Inf, -Inf)) {
    g <- gram
    g[1, 2] <- bad
    expect_error(fit(g), "Gram matrix must be finite")
    expect_identical(fit(g, list(3L)), fit(gram, list(3L)))
    g <- gram
    g[1, 1] <- bad
    expect_error(fit(g), "cross-products must be finite")
    cross <- xty
    cross[2, 1] <- bad
    expect_error(fit(cross = cross), "cross-products must be finite")
  }
})


test_that("shared Gram scans preserve unused nonfinite values and late errors", {
  set.seed(10141)
  n <- 60L
  p <- 33L
  x <- scale(matrix(rnorm(n * p), n))
  y <- matrix(rnorm(n * 12L), n)
  gram <- crossprod(x)
  xty <- crossprod(x, y)
  candidates <- rep(list(32:1), ncol(y))
  run <- function(g, c = candidates, z = xty) {
    fit_greedy_l0_batch(g, z, colSums(y^2), c, n)
  }
  reference <- run(gram)
  for (bad in c(NA_real_, NaN, Inf, -Inf)) {
    g <- gram
    g[33L, 33L] <- bad
    expect_identical(run(g), reference)
    g <- gram
    g[33L, 1L] <- bad
    expect_identical(run(g), reference)
    late <- candidates
    late[[12L]] <- 1:33
    expect_error(run(g, late), "Gram matrix must be finite")
  }
  for (bad in list(c(1L, 1L), c(1L, NA_integer_), 0L, -1L, 34L)) {
    late <- candidates
    late[[12L]] <- bad
    expect_error(run(gram, late), "unique and one-based")
  }
  cross <- xty
  cross[1L, 12L] <- Inf
  expect_error(run(gram, z = cross), "cross-products must be finite")
  for (i in seq_len(ncol(y))) {
    single <- fit_greedy_l0_batch(gram, xty[, i, drop = FALSE], colSums(y^2)[i], candidates[i], n)
    selected <- reference$target_index == i
    for (field in c("predictor_index", "standardized_beta", "deletion_delta_bic")) {
      expect_identical(reference[[field]][selected], single[[field]])
    }
    for (field in c("candidate_size", "support_size", "rss", "bic", "status")) {
      expect_identical(reference[[field]][i], single[[field]])
    }
  }
})
