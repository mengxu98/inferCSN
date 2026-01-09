# Construct network for single target gene

Construct network for single target gene

## Usage

``` r
single_network(
  matrix,
  regulators,
  target,
  cross_validation = FALSE,
  seed = 1,
  penalty = "L0",
  r_squared_threshold = 0,
  n_folds = 5,
  verbose = TRUE,
  ...
)
```

## Arguments

- matrix:

  An expression matrix.

- regulators:

  The regulator genes for which to infer the regulatory network.

- target:

  The target gene.

- cross_validation:

  Whether to use cross-validation. Default is `FALSE`.

- seed:

  The random seed for cross-validation. Default is `1`.

- penalty:

  The type of regularization, default is `"L0"`. This can take either
  one of the following choices: `"L0"`, `"L0L1"`, and `"L0L2"`. For
  high-dimensional and sparse data, `"L0L2"` is more effective.

- r_squared_threshold:

  Threshold of \\R^2\\ coefficient. Default is `0`.

- n_folds:

  The number of folds for cross-validation. Default is `5`.

- verbose:

  Whether to print progress messages. Default is `TRUE`.

- ...:

  Parameters for other methods.

## Value

A data frame of the single target gene network. The data frame has three
columns: regulator, target, and weight.

## Examples

``` r
data(example_matrix)
head(
  single_network(
    example_matrix,
    regulators = colnames(example_matrix),
    target = "g1"
  )
)
#>   regulator target       weight
#> 1       g10     g1  0.009932787
#> 2       g11     g1  0.010286958
#> 3       g12     g1  0.017584462
#> 4       g13     g1 -0.002483363
#> 5       g14     g1 -0.021346832
#> 6       g15     g1  0.027390264
head(
  single_network(
    example_matrix,
    regulators = colnames(example_matrix),
    target = "g1",
    cross_validation = TRUE
  )
)
#>   regulator target    weight
#> 1       g10     g1 0.0000000
#> 2       g11     g1 0.0000000
#> 3       g12     g1 0.0220384
#> 4       g13     g1 0.0000000
#> 5       g14     g1 0.0000000
#> 6       g15     g1 0.0000000

single_network(
  example_matrix,
  regulators = c("g1", "g2", "g3"),
  target = "g1"
)
#>   regulator target     weight
#> 1        g2     g1  0.9876086
#> 2        g3     g1 -0.1569370
single_network(
  example_matrix,
  regulators = c("g1", "g2"),
  target = "g1"
)
#> ! [2026-01-09 07:09:54] Less than 2 regulators found when modeling: "g1"
#> NULL
```
