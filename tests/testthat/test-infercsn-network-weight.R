test_that("inferCSN exports only target-specific greedy-L0 support", {
  set.seed(2026)
  n <- 160
  tf_a <- rnorm(n)
  tf_b <- rnorm(n)
  expression <- cbind(
    tf_a = tf_a,
    tf_b = tf_b,
    target_a = 2 * tf_a + rnorm(n, sd = 0.03),
    target_b = -1.4 * tf_b + rnorm(n, sd = 0.03),
    noise = rnorm(n)
  )

  network <- inferCSN(
    expression,
    regulators = c("tf_a", "tf_b"),
    targets = c("target_a", "target_b"),
    cores = 1,
    verbose = FALSE
  )

  expect_named(network, c("regulator", "target", "weight"))
  expect_setequal(
    paste(network$regulator, network$target, sep = "->"),
    c("tf_a->target_a", "tf_b->target_b")
  )
})

test_that("inferCSN defaults to no user-imposed support cap", {
  set.seed(2041)
  n <- 300L
  p <- 24L
  regulators <- matrix(
    rnorm(n * p),
    nrow = n,
    dimnames = list(NULL, paste0("tf_", seq_len(p)))
  )
  target <- as.numeric(
    regulators %*% seq(0.8, 1.3, length.out = p) + rnorm(n, sd = 0.02)
  )
  expression <- cbind(regulators, target = target)

  network <- inferCSN(
    expression,
    regulators = colnames(regulators),
    targets = "target",
    cores = 1,
    verbose = FALSE
  )
  cap20_network <- inferCSN(
    expression,
    regulators = colnames(regulators),
    targets = "target",
    max_support_size = 20,
    cores = 1,
    verbose = FALSE
  )
  capped_network <- inferCSN(
    expression,
    regulators = colnames(regulators),
    targets = "target",
    max_support_size = 5,
    cores = 1,
    verbose = FALSE
  )

  expect_gt(nrow(network), 20L)
  expect_lte(nrow(cap20_network), 20L)
  expect_gt(nrow(network), nrow(cap20_network))
  expect_lte(nrow(capped_network), 5L)
})

test_that("inferCSN weights are invariant to positive gene rescaling", {
  set.seed(2027)
  n <- 180
  tf_a <- rnorm(n)
  tf_b <- rnorm(n)
  expression <- cbind(
    tf_a = tf_a,
    tf_b = tf_b,
    target_a = 2 * tf_a + rnorm(n, sd = 0.02),
    target_b = -1.4 * tf_b + rnorm(n, sd = 0.02),
    noise = rnorm(n)
  )
  scaled_expression <- expression
  scaled_expression[, "tf_a"] <- 100 * scaled_expression[, "tf_a"]

  infer <- function(x) {
    inferCSN(
      x,
      regulators = c("tf_a", "tf_b"),
      targets = c("target_a", "target_b"),
      cores = 1,
      verbose = FALSE
    )
  }
  network <- infer(expression)
  scaled_network <- infer(scaled_expression)

  expect_equal(scaled_network, network, tolerance = 1e-12)
})

test_that("pseudotime-lag weights use the same scale-invariant network ranking", {
  set.seed(2028)
  n <- 200
  lag <- 10
  tf_a <- rnorm(n)
  tf_b <- rnorm(n)
  target_a <- rnorm(n, sd = 0.1)
  target_b <- rnorm(n, sd = 0.1)
  target_a[(lag + 1):n] <-
    1.8 * tf_a[seq_len(n - lag)] + rnorm(n - lag, sd = 0.03)
  target_b[(lag + 1):n] <-
    -1.2 * tf_b[seq_len(n - lag)] + rnorm(n - lag, sd = 0.06)
  expression <- cbind(
    tf_a = tf_a,
    tf_b = tf_b,
    target_a = target_a,
    target_b = target_b,
    noise = rnorm(n)
  )
  scaled_expression <- expression
  scaled_expression[, "tf_a"] <- 100 * scaled_expression[, "tf_a"]
  scaled_expression[, "target_b"] <- 0.01 * scaled_expression[, "target_b"]

  infer <- function(x) {
    inferCSN(
      x,
      pseudotime = seq_len(n),
      regulators = c("tf_a", "tf_b"),
      targets = c("target_a", "target_b"),
      cores = 1,
      verbose = FALSE
    )
  }
  network <- infer(expression)
  scaled_network <- infer(scaled_expression)

  expect_equal(scaled_network, network, tolerance = 1e-12)
  expect_equal(network$weight, c(0.75, -0.25), tolerance = 1e-12)
})

test_that("network weights and ordering are deterministic across core counts", {
  set.seed(2029)
  n <- 220
  regulators <- replicate(4, rnorm(n))
  colnames(regulators) <- paste0("tf_", seq_len(ncol(regulators)))
  expression <- cbind(
    regulators,
    target_a = 1.5 * regulators[, 1] - 0.8 * regulators[, 2] + rnorm(n, sd = 0.1),
    target_b = -1.1 * regulators[, 3] + rnorm(n, sd = 0.2),
    target_c = 0.7 * regulators[, 4] + rnorm(n, sd = 0.3)
  )
  infer <- function(cores) {
    inferCSN(
      expression,
      regulators = colnames(regulators),
      targets = c("target_a", "target_b", "target_c"),
      cores = cores,
      verbose = FALSE
    )
  }

  single_core <- infer(1)
  two_cores <- infer(2)

  expect_equal(two_cores, single_core, tolerance = 1e-12)
  expect_true(all(diff(abs(single_core$weight)) <= 0))
})

test_that("single_network is the same support-only inference for one target", {
  set.seed(2034)
  n <- 140L
  tf_a <- rnorm(n)
  tf_b <- rnorm(n)
  expression <- cbind(
    tf_a = tf_a,
    tf_b = tf_b,
    target = 1.7 * tf_a - 0.6 * tf_b + rnorm(n, sd = 0.05)
  )

  direct <- inferCSN(
    expression,
    regulators = c("tf_a", "tf_b"),
    targets = "target",
    cores = 1,
    verbose = FALSE
  )
  delegated <- single_network(
    expression,
    regulators = c("tf_a", "tf_b"),
    target = "target",
    cores = 1,
    verbose = FALSE
  )

  expect_equal(delegated, direct, tolerance = 1e-12)
})
