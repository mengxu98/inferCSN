test_that("candidate limits use numeric significance, with stable ties", {
  statistics <- data.frame(padjust = c(.04, .001, .001, NA, .5), row.names = c("a", "b", "c", "d", "e"))
  expect_identical(select_trend_features(statistics = statistics)$features, c("a", "b", "c"))
  expect_identical(select_trend_features(statistics = statistics, n_candidates = 2)$features, c("b", "c"))
  expect_length(select_trend_features(statistics = statistics, padjust_threshold = 0)$features, 0)
  expect_error(select_trend_features(statistics = statistics, n_candidates = 0), "positive")
  expect_error(select_trend_features(statistics = statistics, cores = 2), "precomputed")
})

test_that("network candidate fitting uses observation by feature input", {
  set.seed(1)
  time <- seq(0, 1, length.out = 80)
  x <- cbind(trend = sin(time * 4) + rnorm(80, sd = .1), noise = rnorm(80))
  out <- select_trend_features(x, time)
  expect_true("trend" %in% out$features)
  expect_equal(out$statistics, thisutils::fit_trends(t(x), time, method = "pretsa")$statistics)
})
