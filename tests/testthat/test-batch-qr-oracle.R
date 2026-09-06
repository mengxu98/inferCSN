
test_that("batch supports and deletion evidence agree with direct QR refits", {
  for (seed in 3101:3108) {
    set.seed(seed)
    n <- 70L
    x <- matrix(rnorm(n * 7), nrow = n)
    x[, 6] <- x[, 1] + x[, 2] + rnorm(n, sd = 0.15)
    x[, 7] <- x[, 3] + rnorm(n, sd = 0.02)
    x <- scale(x)
    y <- scale(cbind(1.2*x[, 1] - 0.7*x[, 2] + rnorm(n, sd = 0.3),
                     x[, 3] + 0.5*x[, 4] + rnorm(n, sd = 0.2)))
    candidates <- list(seq_len(7), c(1L, 3L, 4L, 5L, 7L))
    cap <- 4L
    fit <- fit_greedy_l0_batch(crossprod(x), crossprod(x, y), colSums(y^2),
                               candidates, n_obs = n, max_support_size = cap)
    for (target in 1:2) {
      keep <- which(fit$target_index == target)
      selected <- fit$predictor_index[keep]
      qr_fit <- function(s) {
        if (!length(s)) return(list(rss = sum(y[, target]^2), coefficient = numeric()))
        model <- lm.fit(x[, s, drop = FALSE], y[, target])
        if (model$rank != length(s)) return(list(rss = Inf, coefficient = numeric()))
        list(rss = sum(model$residuals^2), coefficient = model$coefficients)
      }
      bic <- function(s) n * log(max(qr_fit(s)$rss / n, 1e-12)) + length(s)*log(n)
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

test_that("constant, singleton, duplicate and wide predictors obey support contracts", {
  set.seed(3190)
  signal <- rnorm(12)
  fixtures <- list(matrix(1, 12, 2), matrix(signal, ncol = 1),
                   cbind(signal, signal, rnorm(12)), matrix(rnorm(12*20), 12, 20))
  for (x in fixtures) {
    fit <- fit_greedy_l0(x, signal + rnorm(12, sd = 0.1), verbose = FALSE)
    expect_lte(length(fit$model$support), min(ncol(x), nrow(x)-2L))
    expect_true(all(is.finite(fit$coefficients$coefficient)))
    selected <- fit$coefficients$coefficient != 0
    expect_true(all(fit$coefficients$deletion_delta_bic[selected] >= 0))
  }
  empty <- inferCSN(cbind(a=rep(1,20),b=rep(2,20)),verbose=FALSE)
  expect_equal(nrow(empty), 0L)
})

test_that("static network output is invariant to core count and row ordering", {
  set.seed(3191)
  x <- matrix(rnorm(80*5),80,5,dimnames=list(paste0("c",1:80),paste0("g",1:5)))
  x[, 4] <- x[, 1] - x[, 2] + rnorm(80,sd=0.1)
  one <- inferCSN(x, cores=1, verbose=FALSE)
  two <- inferCSN(x, cores=2, verbose=FALSE)
  expect_identical(one, two)
  canonical <- function(z) {
    z <- z[order(z$regulator, z$target), , drop = FALSE]
    rownames(z) <- NULL
    z
  }
  expect_equal(canonical(one), canonical(inferCSN(x[80:1,],cores=1,verbose=FALSE)), tolerance=1e-12)
})

test_that("batch agrees with a QR backward reference removing multiple redundant predictors", {
  set.seed(3192)
  n <- 100L
  x <- scale(matrix(rnorm(n*6),n,6))
  x[, 6] <- scale(x[, 1] + x[, 2] + rnorm(n, sd=.4))
  y <- as.numeric(scale(1.5*x[, 1] - x[, 2] + rnorm(n, sd=.01)))
  bic <- function(s) {
    rss <- if (length(s)) sum(lm.fit(x[,s,drop=FALSE],y)$residuals^2) else sum(y^2)
    n*log(max(rss/n,1e-12)) + length(s)*log(n)
  }
  support <- 1:6
  trajectory <- bic(support)
  repeat {
    trials <- lapply(seq_along(support), function(i) support[-i])
    scores <- vapply(trials,bic,numeric(1))
    if (!length(scores) || min(scores) >= tail(trajectory,1)-1e-8) break
    support <- trials[[which.min(scores)]]
    trajectory <- c(trajectory,min(scores))
  }
  expect_gte(length(trajectory)-1L, 2L)
  expect_true(all(diff(trajectory)<0))
  fit <- fit_greedy_l0_batch(crossprod(x),matrix(crossprod(x,y),ncol=1),sum(y^2),list(1:6),n)
  expect_setequal(fit$predictor_index,support)
  expect_equal(fit$bic,tail(trajectory,1),tolerance=1e-7)
})
