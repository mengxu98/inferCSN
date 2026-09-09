test_that("public static weights reconstruct from complete deletion evidence", {
  set.seed(20260815)
  n <- 80L
  tf1 <- stats::rnorm(n)
  tf2 <- stats::rnorm(n)
  object <- cbind(
    TF1 = tf1,
    TF2 = tf2,
    G1 = 1.4 * tf1 - 0.4 * tf2 + stats::rnorm(n, sd = 0.3),
    G2 = 0.8 * tf2 + stats::rnorm(n, sd = 0.4)
  )
  rownames(object) <- paste0("cell", seq_len(n))
  params <- list(
    min_improvement = 1e-10,
    pseudotime_lag_fraction = 0.05,
    pseudotime_lag_steps = NULL,
    regulators = c("TF1", "TF2"),
    targets = c("G1", "G2"),
    max_support_size = NULL,
    cores = 1L
  )
  archived <- inferCSN:::infer_network(
    t(object), colnames(object), matrix(numeric(0), nrow(object), 0L), params
  )
  public <- inferCSN(
    object,
    regulators = params$regulators, targets = params$targets,
    cores = 1L, verbose = FALSE
  )
  group_rank <- rank(-archived$deletion_delta_bic, ties.method = "min")
  group_index <- match(group_rank, sort(unique(group_rank)))
  expected <- sign(archived$standardized_beta) *
    (1 - (group_index - 0.5) / length(unique(group_index)))
  expect_identical(names(public), c("regulator", "target", "weight"))
  expect_equal(archived$weight, expected, tolerance = 0)
  expect_identical(public, archived[, c("regulator", "target", "weight")])
})

test_that("batch Gram interface matches the single-response solver", {
  set.seed(20260830)
  x <- scale(matrix(stats::rnorm(240), nrow = 80L, ncol = 3L))
  y <- scale(cbind(
    1.2 * x[, 1L] + stats::rnorm(80L, sd = 0.2),
    -0.8 * x[, 2L] + stats::rnorm(80L, sd = 0.2)
  ))
  batch <- fit_greedy_l0_batch(
    crossprod(x), crossprod(x, y), colSums(y^2),
    list(1:3, 1:3),
    n_obs = nrow(x)
  )
  singles <- lapply(seq_len(ncol(y)), function(index) {
    fit_greedy_l0(x, y[, index], verbose = FALSE)
  })
  for (index in seq_along(singles)) {
    selected <- singles[[index]]$coefficients$coefficient != 0
    observed <- batch$target_index == index
    expect_equal(
      batch$standardized_beta[observed],
      unname(singles[[index]]$coefficients$coefficient[selected]),
      tolerance = 1e-10
    )
    expect_equal(
      batch$deletion_delta_bic[observed],
      singles[[index]]$coefficients$deletion_delta_bic[selected],
      tolerance = 1e-8
    )
  }
})

test_that("lagged deletion-evidence weights are deterministic across cores", {
  set.seed(20260816)
  n <- 60L
  tf1 <- stats::rnorm(n)
  tf2 <- stats::rnorm(n)
  object <- cbind(
    TF1 = tf1,
    TF2 = tf2,
    G1 = c(0, head(tf1, -1)) + stats::rnorm(n, sd = 0.2),
    G2 = c(0, head(tf2, -1)) + stats::rnorm(n, sd = 0.2)
  )
  rownames(object) <- paste0("cell", seq_len(n))
  one <- inferCSN(
    object,
    pseudotime = seq_len(n), lag_steps = 1L,
    regulators = c("TF1", "TF2"), targets = c("G1", "G2"),
    cores = 1L, verbose = FALSE
  )
  two <- inferCSN(
    object,
    pseudotime = seq_len(n), lag_steps = 1L,
    regulators = c("TF1", "TF2"), targets = c("G1", "G2"),
    cores = 2L, verbose = FALSE
  )
  expect_identical(one, two)
})

test_that("fit_greedy_l0 archives deletion evidence for selected support", {
  set.seed(20260817)
  x <- cbind(a = stats::rnorm(70), b = stats::rnorm(70))
  y <- 1.2 * x[, "a"] + stats::rnorm(70, sd = 0.25)
  fit <- fit_greedy_l0(x, y, verbose = FALSE)
  selected <- fit$coefficients$coefficient != 0
  expect_true(all(is.finite(fit$coefficients$deletion_delta_bic[selected])))
  expect_true(all(fit$coefficients$deletion_delta_bic[selected] >= 0))
  expect_true(all(is.na(fit$coefficients$deletion_delta_bic[!selected])))
})

test_that("batch solver preserves analytic SSE near a perfect fit", {
  gram <- matrix(c(
    799.000000000005, 658.642698426331, 104.012480436715,
    -2.47176430059963, -1.3911468059158, -2.75562303346975,
    -1.3911468059158, 176.879813218877, -1.39114680591578,
    658.642698426331, 799.000000000003, -4.0595572761924,
    -1.77678178182818, -1, -1.98082835093435, -1.00000000000001,
    -1.60310235263307, -1, 104.012480436715, -4.0595572761924,
    799.000000000002, 114.228237804294, 89.851188663258,
    49.8353973720042, 89.851188663258, 114.663194811407,
    188.731680519065, -2.47176430059963, -1.77678178182818,
    114.228237804294, 799.000000000009, 276.71259097775,
    119.82363936541, 276.71259097775, 672.268756438507,
    134.156368828086, -1.3911468059158, -1, 89.851188663258,
    276.71259097775, 799.000000000004, -1.98082835093434,
    799.000000000004, -1.60310235263307, -1, -2.75562303346975,
    -1.98082835093435, 49.8353973720042, 119.82363936541,
    -1.98082835093434, 798.999999999998, -1.98082835093434,
    83.0007074957214, 86.6657262156137, -1.3911468059158,
    -1.00000000000001, 89.851188663258, 276.71259097775,
    799.000000000004, -1.98082835093434, 799.000000000004,
    -1.60310235263307, -1, 176.879813218877, -1.60310235263307,
    114.663194811407, 672.268756438507, -1.60310235263307,
    83.0007074957214, -1.60310235263307, 799.000000000002,
    -1.60310235263307, -1.39114680591578, -1, 188.731680519065,
    134.156368828086, -1, 86.6657262156137, -1,
    -1.60310235263307, 799.000000000001
  ), nrow = 9L)
  xty <- matrix(c(
    -1.39114680591578, -1, 118.84240286243, 362.685127742302,
    -1, 235.696974341759, -1, 288.457398405829, -1
  ), ncol = 1L)

  fit <- fit_greedy_l0_batch(
    gram, xty,
    response_ss = 799, candidates = list(1:9), n_obs = 800
  )

  expect_equal(fit$support_size, 6L)
  expect_true(all(fit$deletion_delta_bic >= 0))
})
