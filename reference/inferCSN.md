# inferring cell-type specific gene regulatory network

inferring cell-type specific gene regulatory network

## Usage

``` r
inferCSN(
  object,
  pseudotime = NULL,
  regulators = NULL,
  targets = NULL,
  max_support_size = NULL,
  lag_fraction = 0.05,
  lag_steps = NULL,
  cores = 1,
  verbose = TRUE,
  method = c("greedy_l0", "L0", "L0L1", "L0L2"),
  ...
)

# S4 method for class 'matrix'
inferCSN(
  object,
  pseudotime = NULL,
  regulators = NULL,
  targets = NULL,
  max_support_size = NULL,
  lag_fraction = 0.05,
  lag_steps = NULL,
  cores = 1,
  verbose = TRUE,
  method = c("greedy_l0", "L0", "L0L1", "L0L2"),
  ...
)

# S4 method for class 'sparseMatrix'
inferCSN(
  object,
  pseudotime = NULL,
  regulators = NULL,
  targets = NULL,
  max_support_size = NULL,
  lag_fraction = 0.05,
  lag_steps = NULL,
  cores = 1,
  verbose = TRUE,
  method = c("greedy_l0", "L0", "L0L1", "L0L2"),
  ...
)
```

## Arguments

- object:

  Numeric expression matrix with cells in rows and genes in columns.

- pseudotime:

  Optional pseudotime vector or branch matrix for either method.
  Regulators in earlier states predict targets in later states; tied
  states are averaged and shared branch transitions are counted once.

- regulators, targets:

  Optional gene subsets.

- max_support_size:

  Optional support-size limit.

- lag_fraction:

  Fractional lag used when `lag_steps` is `NULL`.

- lag_steps:

  Optional integer lag. L0Learn centers and scales lagged expression
  within each branch and defaults to fitting without an intercept.

- cores:

  Number of inference workers.

- verbose:

  Whether to report progress.

- method:

  `greedy_l0` (default), or L0Learn with the `L0`, `L0L1`, or `L0L2`
  penalty.

- ...:

  Arguments passed to the method.

## Value

A data frame containing exactly `regulator`, `target`, and `weight`.

## Examples

``` r
data(example_matrix)
data(example_meta_data)
network_table <- inferCSN(
  example_matrix,
  pseudotime = example_meta_data$pseudotime
)
#> ℹ [2026-10-07 03:30:41] Inferring network for <matrix/array>...
#> ✔ [2026-10-07 03:30:41] Inferring network done
#> ℹ [2026-10-07 03:30:41] Network information:
#> ℹ                         Edges Regulators Targets
#> ℹ                       1    11          6       6
head(network_table)
#>   regulator target     weight
#> 1        g5     g1 -0.9545455
#> 2        g4     g5  0.8636364
#> 3        g3     g4  0.7727273
#> 4        g1     g2 -0.6818182
#> 5        g1     g6  0.5909091
#> 6        g6     g3  0.5000000

inferCSN(
  example_matrix,
  regulators = c("g1", "g2"),
  targets = c("g3", "g4")
)
#> ℹ [2026-10-07 03:30:41] Inferring network for <matrix/array>...
#> ✔ [2026-10-07 03:30:41] Inferring network done
#> ℹ [2026-10-07 03:30:41] Network information:
#> ℹ                         Edges Regulators Targets
#> ℹ                       1     4          2       2
#>   regulator target weight
#> 1        g2     g3 -0.875
#> 2        g1     g4 -0.625
#> 3        g2     g4 -0.375
#> 4        g1     g3 -0.125
```
