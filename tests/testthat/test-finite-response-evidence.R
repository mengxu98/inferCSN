test_that("deletion evidence uses the same finite-response rows as the fit", {
  set.seed(1048)
  x <- cbind(a = rnorm(60), b = rnorm(60), c = rnorm(60))
  y <- 1.4 * x[, 1] - 0.6 * x[, 2] + rnorm(60, sd = 0.25)
  y[c(2, 29)] <- NA_real_
  x[c(2, 29), 1] <- c(1000, -1000)
  keep <- is.finite(y)
  observed <- fit_greedy_l0(x, y, verbose = FALSE)
  expected <- fit_greedy_l0(x[keep, , drop = FALSE], y[keep], verbose = FALSE)
  expect_identical(observed$model, expected$model)
  expect_identical(observed$coefficients$coefficient, expected$coefficients$coefficient)
  expect_equal(observed$coefficients$deletion_delta_bic,
               expected$coefficients$deletion_delta_bic, tolerance = 1e-12)
})
