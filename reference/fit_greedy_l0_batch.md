# Fit batched greedy L0 regressions

Fits multiple responses from shared predictor cross-products.

## Usage

``` r
fit_greedy_l0_batch(
  gram,
  xty,
  response_ss,
  candidates,
  n_obs,
  max_support_size = NULL,
  min_improvement = 1e-10
)
```

## Arguments

- gram:

  Numeric square predictor cross-product matrix.

- xty:

  Numeric predictor-by-response cross-product matrix.

- response_ss:

  Numeric response sums of squares.

- candidates:

  List of one-based predictor indices, one element per response.

- n_obs:

  Number of observations used for the cross-products.

- max_support_size:

  Optional positive integer support limit.

- min_improvement:

  Minimum RSS reduction proposed during forward search.

## Value

Batched coefficients, deletion evidence, and fit summaries.
