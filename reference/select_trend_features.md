# Select features with significant trends

Select features with significant trends

## Usage

``` r
select_trend_features(
  x = NULL,
  pseudotime = NULL,
  method = c("pretsa", "gam"),
  padjust_threshold = 0.05,
  n_candidates = NULL,
  statistics = NULL,
  ...
)
```

## Arguments

- x:

  Numeric matrix with observations in rows and features in columns, or
  NULL when supplying `statistics`.

- pseudotime:

  Numeric coordinate per observation.

- method:

  Fitting method: `"pretsa"` or `"gam"`.

- padjust_threshold:

  Strict adjusted P-value cutoff between zero and one.

- n_candidates:

  Maximum number of selected features, or NULL for no cap.

- statistics:

  Precomputed table with feature row names and `padjust` (or `pvalue` if
  `padjust` is absent).

- ...:

  Numerical options passed to
  [`thisutils::fit_trends()`](https://mengxu98.github.io/thisutils/reference/fit_trends.html).

## Value

A list with `features`, `statistics`, and `fit` (NULL for precomputed
input). Nonfinite P-values are excluded. Input order is retained unless
the cap is exceeded; then the smallest P-values are selected with stable
ties.

## Examples

``` r
t <- seq(0, 1, length.out = 40)
x <- cbind(a = sin(t * 5), b = cos(t * 3))
select_trend_features(x, t)$features
#> [1] "a" "b"
```
