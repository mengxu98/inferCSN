# Sparse regression model

Sparse regression model

## Usage

``` r
fit_srm(
  x,
  y,
  cross_validation = FALSE,
  seed = 1,
  penalty = "L0",
  regulators_num = ncol(x),
  n_folds = 5,
  verbose = TRUE,
  ...
)
```

## Arguments

- x:

  The matrix of regulators.

- y:

  The vector of target.

- cross_validation:

  Whether to use cross-validation. Default is `FALSE`.

- seed:

  The random seed for cross-validation. Default is `1`.

- penalty:

  The type of regularization. This can take either one of the following
  choices: `"L0"`, `"L0L1"`, and `"L0L2"`. For high-dimensional and
  sparse data, `"L0L2"` is more effective. Default is `"L0"`.

- regulators_num:

  The number of regulators for target.

- n_folds:

  The number of folds for cross-validation. Default is `5`.

- verbose:

  Whether to print progress messages. Default is `TRUE`.

- ...:

  Parameters for other methods.

## Value

A list of the sparse regression model. The list has the three
components: model, metrics, and coefficients.

## Examples

``` r
data(example_matrix)
fit_srm(
  x = example_matrix[, -1],
  y = example_matrix[, 1]
)
#> $model
#> 
#> $metrics
#> $metrics$r_squared
#> [1] 0.4358336
#> 
#> 
#> $coefficients
#> $coefficients$variable
#>  [1] "g10" "g11" "g12" "g13" "g14" "g15" "g16" "g17" "g18" "g2"  "g3"  "g4" 
#> [13] "g5"  "g6"  "g7"  "g8"  "g9" 
#> 
#> $coefficients$coefficient
#>  [1]  0.006073601  0.006290166  0.010752371 -0.001518502 -0.013052947
#>  [6]  0.016748325 -0.024070066  0.069284763 -0.563969602  0.217868939
#> [11]  0.007553373 -0.017714243  0.023644447 -0.018219713  0.015716106
#> [16]  0.005408364 -0.028352269
#> 
#> 
```
