test_that("L0Learn single-target results contain only selected edges", {
  skip_if_not_installed("L0Learn")
  set.seed(7)
  x <- matrix(rnorm(100 * 4), 100, 4,
    dimnames = list(NULL, c("a", "b", "c", "target"))
  )
  x[, "target"] <- 3 * x[, "a"] + rnorm(100, sd = 0.01)
  result <- single_network(x,
    regulators = c("a", "b", "c"), target = "target",
    method = "L0", max_support_size = 1, verbose = FALSE
  )
  expect_equal(nrow(result), 1L)
  expect_identical(result$regulator, "a")
  expect_true(all(is.finite(result$weight) & result$weight != 0))

  empty <- single_network(x,
    regulators = "a", target = "target", method = "L0", verbose = FALSE
  )
  expect_s3_class(empty, "data.frame")
  expect_identical(names(empty), c("regulator", "target", "weight"))
  expect_equal(nrow(empty), 0L)
})
