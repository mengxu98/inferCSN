# Construct network for single target gene

Construct network for single target gene

## Usage

``` r
single_network(
  matrix,
  regulators,
  target,
  pseudotime = NULL,
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

- matrix:

  An expression matrix.

- regulators:

  Candidate regulator genes.

- target:

  The target gene.

- pseudotime:

  Optional pseudotime vector or branch matrix passed to \[inferCSN()\].

- max_support_size:

  Optional support-size cap passed to \[inferCSN()\].

- lag_fraction:

  Fractional state lag passed to \[inferCSN()\].

- lag_steps:

  Optional integer state lag passed to \[inferCSN()\].

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

A data frame with regulator, target, and weight columns. Greedy-L0
returns selected edges. L0Learn retains the original per-regulator
coefficients, including zeros; \[inferCSN()\] removes zero weights from
the complete network.

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
#> ℹ [2026-10-08 08:44:54] Inferring network for <matrix/array>...
#> ✔ [2026-10-08 08:44:54] Inferring network done
#> ℹ [2026-10-08 08:44:54] Network information:
#> ℹ                         Edges Regulators Targets
#> ℹ                       1     2          2       1
#>   regulator target weight
#> 1        g6     g1   0.75
#> 2        g5     g1  -0.25
single_network(
  example_matrix,
  regulators = c("g1", "g2", "g3"),
  target = "g1"
)
#> ℹ [2026-10-08 08:44:54] Inferring network for <matrix/array>...
#> ✔ [2026-10-08 08:44:54] Inferring network done
#> ℹ [2026-10-08 08:44:54] Network information:
#> ℹ                         Edges Regulators Targets
#> ℹ                       1     2          2       1
#>   regulator target weight
#> 1        g2     g1  -0.75
#> 2        g3     g1  -0.25
```
