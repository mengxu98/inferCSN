# Fit greedy L0 regression

Fits a standardized sparse linear model by BIC-decreasing support
updates.

## Usage

``` r
fit_greedy_l0(
  x,
  y,
  max_support_size = NULL,
  min_improvement = 1e-10,
  verbose = TRUE
)
```

## Arguments

- x:

  Numeric predictor matrix with observations in rows.

- y:

  Numeric response vector.

- max_support_size:

  Optional positive integer support limit.

- min_improvement:

  Minimum RSS improvement considered.

- verbose:

  Logical verbosity flag.

## Value

A list with \`model\`, \`metrics\`, and \`coefficients\`.
