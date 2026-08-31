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
  verbose = TRUE
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

## Value

A data frame containing only selected edges for the requested target.
The data frame has three columns: regulator, target, and weight.

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
#> ℹ [2026-08-31 02:40:03] Inferring network for <matrix/array>...
#> ◌ [2026-08-31 02:40:03] Checking parameters...
#> ✔ [2026-08-31 02:40:03] Inferring network done
#> ℹ [2026-08-31 02:40:03] Network information:
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
#> ℹ [2026-08-31 02:40:03] Inferring network for <matrix/array>...
#> ◌ [2026-08-31 02:40:03] Checking parameters...
#> ✔ [2026-08-31 02:40:03] Inferring network done
#> ℹ [2026-08-31 02:40:03] Network information:
#> ℹ                         Edges Regulators Targets
#> ℹ                       1     2          2       1
#>   regulator target weight
#> 1        g2     g1  -0.75
#> 2        g3     g1  -0.25
```
