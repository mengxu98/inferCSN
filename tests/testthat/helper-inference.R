support_neighbours <- function(selected, candidates, cap) {
  outside <- setdiff(candidates, selected)
  neighbours <- lapply(seq_along(selected), function(i) selected[-i])
  if (length(selected) < cap) {
    neighbours <- c(neighbours, lapply(outside, function(j) sort(c(selected, j))))
  }
  for (i in seq_along(selected)) {
    neighbours <- c(neighbours, lapply(outside, function(j) sort(c(selected[-i], j))))
  }
  neighbours
}

expect_network_across_cores <- function(object, ...) {
  one <- inferCSN(object, cores = 1L, verbose = FALSE, ...)
  expect_identical(inferCSN(object, cores = 2L, verbose = FALSE, ...), one)
  one
}
