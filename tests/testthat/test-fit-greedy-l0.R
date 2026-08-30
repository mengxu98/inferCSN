test_that("fit_greedy_l0 recovers a sparse signal and reports its path model", {
  set.seed(11)
  x <- cbind(
    tf_a = rnorm(100),
    tf_b = rnorm(100),
    noise = rnorm(100)
  )
  y <- 1.8 * x[, "tf_a"] - 1.2 * x[, "tf_b"] + rnorm(100, sd = 0.05)

  fit <- fit_greedy_l0(x, y)

  expect_identical(fit$model$algorithm, "greedy_l0")
  expect_setequal(fit$model$support, c("tf_a", "tf_b"))
  expect_gt(fit$metrics$r_squared, 0.95)
  expect_equal(
    unname(fit$coefficients$coefficient[match("noise", fit$coefficients$variable)]),
    0
  )
})

test_that("fit_greedy_l0 retains the empty model when no candidate improves RSS", {
  x <- cbind(
    tf_a = rep(1, 20),
    tf_b = rep(2, 20)
  )
  y <- rep(c(-1, 1), 10)

  fit <- fit_greedy_l0(x, y)

  expect_length(fit$model$support, 0)
  expect_true(all(fit$coefficients$coefficient == 0))
  expect_equal(fit$metrics$r_squared, 0)
})

test_that("fit_greedy_l0 caps support by the number of observations", {
  set.seed(12)
  x <- matrix(rnorm(6 * 12), nrow = 6)
  colnames(x) <- paste0("tf_", seq_len(ncol(x)))
  y <- rowSums(x[, seq_len(6), drop = FALSE])

  fit <- fit_greedy_l0(x, y, max_support_size = 20)

  expect_lte(length(fit$model$support), nrow(x) - 2)
})

test_that("fit_greedy_l0 has no fixed default support cap", {
  expect_null(formals(fit_greedy_l0)$max_support_size)

  set.seed(2040)
  n <- 300L
  p <- 24L
  x <- matrix(
    rnorm(n * p),
    nrow = n,
    dimnames = list(NULL, paste0("tf_", seq_len(p)))
  )
  y <- as.numeric(
    x %*% seq(0.8, 1.3, length.out = p) + rnorm(n, sd = 0.02)
  )

  fit <- fit_greedy_l0(x, y, verbose = FALSE)
  cap20_fit <- fit_greedy_l0(
    x,
    y,
    max_support_size = 20,
    verbose = FALSE
  )
  capped_fit <- fit_greedy_l0(
    x,
    y,
    max_support_size = 5,
    verbose = FALSE
  )

  expect_gt(length(fit$model$support), 20L)
  expect_lte(length(cap20_fit$model$support), 20L)
  expect_gt(length(fit$model$support), length(cap20_fit$model$support))
  expect_lte(length(capped_fit$model$support), 5L)
})

test_that("fit_greedy_l0 returns scale-invariant standardized coefficients", {
  set.seed(13)
  x <- cbind(
    tf_a = rnorm(120),
    tf_b = rnorm(120),
    noise = rnorm(120)
  )
  y <- 1.8 * x[, "tf_a"] - 1.2 * x[, "tf_b"] + rnorm(120, sd = 0.05)
  scaled_x <- x
  scaled_x[, "tf_a"] <- 100 * scaled_x[, "tf_a"]

  fit <- fit_greedy_l0(x, y)
  scaled_fit <- fit_greedy_l0(scaled_x, 0.01 * y)

  expect_equal(
    scaled_fit$coefficients$coefficient,
    fit$coefficients$coefficient,
    tolerance = 1e-12
  )
  expect_equal(scaled_fit$metrics, fit$metrics, tolerance = 1e-12)
})

test_that("fit_greedy_l0 escapes a correlated proxy with a one-variable swap", {
  set.seed(1)
  n <- 160L
  x_1 <- rnorm(n)
  x_2 <- rnorm(n)
  x <- cbind(
    x_1 = x_1,
    x_2 = x_2,
    proxy = x_1 + x_2 + rnorm(n, sd = 0.25)
  )
  y <- x_1 + x_2 + rnorm(n, sd = 0.2)

  fit <- fit_greedy_l0(x, y, max_support_size = 2)

  expect_setequal(fit$model$support, c("x_1", "x_2"))
  expect_gt(fit$metrics$r_squared, 0.95)
})

test_that("swap admissibility is evaluated after deleting the incumbent", {
  set.seed(2)
  n <- 80L
  x_1 <- rnorm(n)
  x_2 <- rnorm(n)
  z <- rnorm(n)
  x <- cbind(
    x_1 = x_1,
    x_2 = x_2,
    proxy = x_1 + x_2 + 1e-5 * z
  )
  y <- x_1 + 2 * x_2 + rnorm(n, sd = 1e-8)

  fit <- fit_greedy_l0(x, y, max_support_size = 2)

  expect_setequal(fit$model$support, c("x_1", "x_2"))
  expect_lt(kappa(crossprod(scale(x[, c("x_1", "x_2")]))), 2)
})

test_that("fit_greedy_l0 never accepts a rank-deficient support", {
  set.seed(2042)
  x_1 <- rnorm(120)
  x_2 <- rnorm(120)
  x <- cbind(
    x_1 = x_1,
    x_1_duplicate = x_1,
    x_2 = x_2
  )
  y <- 1.4 * x_1 - 0.8 * x_2 + rnorm(120, sd = 0.04)

  fit <- fit_greedy_l0(x, y, max_support_size = 3)

  expect_true("x_2" %in% fit$model$support)
  expect_lte(
    sum(c("x_1", "x_1_duplicate") %in% fit$model$support),
    1L
  )
  expect_true(all(is.finite(fit$coefficients$coefficient)))
})

test_that("fit_greedy_l0 terminates at a delete-add-swap BIC local optimum", {
  set.seed(2043)
  n <- 90L
  p <- 8L
  x <- matrix(
    rnorm(n * p),
    nrow = n,
    dimnames = list(NULL, paste0("tf_", seq_len(p)))
  )
  x[, "tf_8"] <- x[, "tf_1"] + x[, "tf_2"] + rnorm(n, sd = 0.35)
  y <- 1.3 * x[, "tf_1"] - 0.9 * x[, "tf_2"] +
    0.6 * x[, "tf_4"] + rnorm(n, sd = 0.3)
  support_cap <- 4L

  fit <- fit_greedy_l0(x, y, max_support_size = support_cap)
  z_x <- scale(x)
  z_y <- as.numeric(scale(y))
  bic <- function(support) {
    rss <- if (!length(support)) {
      sum(z_y^2)
    } else {
      sum(lm.fit(z_x[, support, drop = FALSE], z_y)$residuals^2)
    }
    n * log(max(rss / n, 1e-12)) + length(support) * log(n)
  }

  selected <- match(fit$model$support, colnames(x))
  outside <- setdiff(seq_len(p), selected)
  neighbours <- list()
  if (length(selected)) {
    neighbours <- c(
      neighbours,
      lapply(seq_along(selected), function(i) selected[-i])
    )
  }
  if (length(selected) < support_cap) {
    neighbours <- c(
      neighbours,
      lapply(outside, function(j) sort(c(selected, j)))
    )
  }
  if (length(selected) && length(outside)) {
    neighbours <- c(
      neighbours,
      unlist(lapply(seq_along(selected), function(i) {
        lapply(outside, function(j) sort(c(selected[-i], j)))
      }), recursive = FALSE)
    )
  }

  expect_equal(fit$model$bic, bic(selected), tolerance = 1e-8)
  neighbour_bic <- vapply(neighbours, bic, numeric(1))
  expect_true(all(neighbour_bic >= fit$model$bic - 1e-8))
})

test_that("near-collinear randomized fits are exhaustive one-exchange optima", {
  implementation_gram_is_invertible <- function(gram) {
    k <- nrow(gram)
    if (!k) {
      return(TRUE)
    }
    work <- gram
    matrix_scale <- max(abs(diag(work)))
    pivot_tolerance <- 1e-10 * max(1, matrix_scale)
    for (column in seq_len(k)) {
      candidates <- column:k
      pivot <- candidates[[which.max(abs(work[candidates, column]))]]
      if (abs(work[pivot, column]) <= pivot_tolerance) {
        return(FALSE)
      }
      if (pivot != column) {
        work[c(column, pivot), ] <- work[c(pivot, column), , drop = FALSE]
      }
      work[column, ] <- work[column, ] / work[column, column]
      for (row in setdiff(seq_len(k), column)) {
        factor <- work[row, column]
        if (factor != 0) {
          work[row, ] <- work[row, ] - factor * work[column, ]
        }
      }
    }
    TRUE
  }

  exhaustive_neighbour_bic <- function(x, y, support, support_cap) {
    z_x <- scale(x)
    z_y <- as.numeric(scale(y))
    p <- ncol(z_x)
    selected <- match(support, colnames(z_x))
    outside <- setdiff(seq_len(p), selected)
    neighbours <- list()
    if (length(selected)) {
      neighbours <- c(
        neighbours,
        lapply(seq_along(selected), function(i) selected[-i])
      )
    }
    if (length(selected) < support_cap) {
      neighbours <- c(
        neighbours,
        lapply(outside, function(j) sort(c(selected, j)))
      )
    }
    if (length(selected) && length(outside)) {
      neighbours <- c(
        neighbours,
        unlist(lapply(seq_along(selected), function(i) {
          lapply(outside, function(j) sort(c(selected[-i], j)))
        }), recursive = FALSE)
      )
    }
    bic <- function(candidate) {
      if (!length(candidate)) {
        rss <- sum(z_y^2)
      } else {
        gram <- crossprod(z_x[, candidate, drop = FALSE])
        if (!implementation_gram_is_invertible(gram)) {
          return(Inf)
        }
        rss <- sum(lm.fit(z_x[, candidate, drop = FALSE], z_y)$residuals^2)
      }
      nrow(z_x) * log(max(rss / nrow(z_x), 1e-12)) +
        length(candidate) * log(nrow(z_x))
    }
    c(current = bic(selected), vapply(neighbours, bic, numeric(1)))
  }

  for (seed in 2101:2112) {
    set.seed(seed)
    n <- 90L
    x_1 <- rnorm(n)
    x_2 <- rnorm(n)
    x <- cbind(
      x_1 = x_1,
      x_2 = x_2,
      proxy = x_1 + x_2 + runif(1, 0.7e-5, 1.3e-5) * rnorm(n),
      noise_1 = rnorm(n),
      noise_2 = rnorm(n)
    )
    y <- x_1 + runif(1, 1.5, 2.5) * x_2 + rnorm(n, sd = 1e-6)
    support_cap <- 2L
    fit <- fit_greedy_l0(x, y, max_support_size = support_cap)
    neighbour_bic <- exhaustive_neighbour_bic(
      x,
      y,
      fit$model$support,
      support_cap
    )

    expect_equal(fit$model$bic, neighbour_bic[["current"]], tolerance = 1e-7)
    expect_true(
      all(neighbour_bic[-1] >= fit$model$bic - 1e-7),
      info = paste("seed", seed)
    )
  }
})
