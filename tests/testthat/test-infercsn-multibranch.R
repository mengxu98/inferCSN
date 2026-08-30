test_that("inferCSN pools all pseudotime branches without column-order dependence", {
  set.seed(2030)
  cells_per_branch <- 120L
  lag <- ceiling(0.05 * cells_per_branch)
  n <- 2L * cells_per_branch

  tf_a <- rnorm(n)
  tf_b <- rnorm(n)
  target_a <- rnorm(n, sd = 0.05)
  target_b <- rnorm(n, sd = 0.05)

  branch_a <- seq_len(cells_per_branch)
  branch_b <- cells_per_branch + seq_len(cells_per_branch)
  target_a[branch_a[(lag + 1L):cells_per_branch]] <-
    2 * tf_a[branch_a[seq_len(cells_per_branch - lag)]] +
    rnorm(cells_per_branch - lag, sd = 0.02)
  target_b[branch_b[(lag + 1L):cells_per_branch]] <-
    -1.6 * tf_b[branch_b[seq_len(cells_per_branch - lag)]] +
    rnorm(cells_per_branch - lag, sd = 0.02)

  expression <- cbind(
    tf_a = tf_a,
    tf_b = tf_b,
    target_a = target_a,
    target_b = target_b
  )
  pseudotime <- cbind(
    branch_a = c(seq_len(cells_per_branch), rep(NA_real_, cells_per_branch)),
    branch_b = c(rep(NA_real_, cells_per_branch), seq_len(cells_per_branch))
  )
  infer <- function(pt) {
    inferCSN(
      expression,
      pseudotime = pt,
      regulators = c("tf_a", "tf_b"),
      targets = c("target_a", "target_b"),
      max_support_size = 1,
      cores = 1,
      verbose = FALSE
    )
  }

  network <- infer(pseudotime)
  reversed_network <- infer(pseudotime[, 2:1, drop = FALSE])

  expect_setequal(
    paste(network$regulator, network$target, sep = "->"),
    c("tf_a->target_a", "tf_b->target_b")
  )
  expect_equal(reversed_network, network, tolerance = 1e-12)
  expect_equal(
    infer(pseudotime[, 1]),
    infer(pseudotime[, 1, drop = FALSE]),
    tolerance = 1e-12
  )
  expect_error(
    infer(pseudotime[-1, , drop = FALSE]),
    "one row.*per cell"
  )
})

test_that("exact duplicate pseudotime branches do not change the network", {
  set.seed(2031)
  n <- 180L
  lag <- ceiling(0.05 * n)
  tf <- rnorm(n)
  target <- rnorm(n, sd = 0.05)
  target[(lag + 1L):n] <-
    1.7 * tf[seq_len(n - lag)] + rnorm(n - lag, sd = 0.02)
  expression <- cbind(tf = tf, target = target, noise = rnorm(n))
  pseudotime <- seq_len(n)

  infer <- function(pt) {
    inferCSN(
      expression,
      pseudotime = pt,
      regulators = c("tf", "noise"),
      targets = "target",
      cores = 1,
      verbose = FALSE
    )
  }

  expect_equal(
    infer(cbind(branch_a = pseudotime, branch_b = pseudotime)),
    infer(pseudotime),
    tolerance = 1e-12
  )
})

test_that("tied pseudotime states are invariant to cell row order", {
  set.seed(2032)
  cells_per_state <- 6L
  n_states <- 20L
  n <- cells_per_state * n_states
  state <- rep(seq_len(n_states), each = cells_per_state)
  tf_state <- rnorm(n_states)
  tf <- rep(tf_state, each = cells_per_state) + rnorm(n, sd = 0.03)
  target_state <- c(0, 1.8 * tf_state[-n_states])
  target <- rep(target_state, each = cells_per_state) + rnorm(n, sd = 0.03)
  expression <- cbind(tf = tf, target = target, noise = rnorm(n))

  infer <- function(x, pt) {
    inferCSN(
      x,
      pseudotime = pt,
      regulators = c("tf", "noise"),
      targets = "target",
      cores = 1,
      verbose = FALSE
    )
  }
  reference <- infer(expression, state)
  permutation <- sample.int(n)
  permuted <- infer(expression[permutation, , drop = FALSE], state[permutation])

  expect_equal(permuted, reference, tolerance = 1e-12)
  expect_identical(
    paste(reference$regulator, reference$target, sep = "->"),
    "tf->target"
  )
})

test_that("overlapping branches use unique transitions for the support bound", {
  set.seed(2033)
  n <- 10L
  p <- 12L
  regulators <- matrix(
    rnorm(n * p),
    nrow = n,
    dimnames = list(NULL, paste0("tf_", seq_len(p)))
  )
  target <- c(0, rowSums(regulators[seq_len(n - 1L), , drop = FALSE]))
  expression <- cbind(regulators, target = target)
  pseudotime <- cbind(
    full = seq_len(n),
    overlapping_subset = c(seq_len(n - 2L), NA_real_, NA_real_)
  )

  network <- inferCSN(
    expression,
    pseudotime = pseudotime,
    regulators = colnames(regulators),
    targets = "target",
    lag_steps = 1,
    cores = 1,
    verbose = FALSE
  )

  unique_transitions <- n - 1L
  distinct_branches <- 2L
  expect_lte(
    nrow(network),
    unique_transitions - distinct_branches - 1L
  )
})

test_that("inferCSN exposes only controls used by the native algorithm", {
  public_controls <- names(formals(methods::getGeneric("inferCSN")))

  expect_setequal(
    public_controls,
    c(
      "object", "pseudotime", "regulators", "targets",
      "max_support_size", "lag_fraction", "lag_steps", "cores", "verbose", "..."
    )
  )
  expect_false(any(c(
    "algorithm", "penalty", "cross_validation", "seed", "n_folds",
    "subsampling_method", "subsampling_ratio"
  ) %in% public_controls))
})
