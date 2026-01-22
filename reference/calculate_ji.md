# Calculate Jaccard Index

Calculates Jaccard Index (JI) metric

## Usage

``` r
calculate_ji(network_table, ground_truth)
```

## Arguments

- network_table:

  A data frame of predicted network structure

- ground_truth:

  A data frame of ground truth network

## Value

A list containing the metric

## Examples

``` r
data(example_matrix)
data("example_ground_truth")
network_table <- inferCSN(example_matrix)
#> ℹ [2026-01-22 03:00:59] Inferring network for <dense matrix>...
#> ◌ [2026-01-22 03:00:59] Checking parameters...
#> ℹ [2026-01-22 03:00:59] Using L0 sparse regression model
#> ℹ [2026-01-22 03:00:59] Using 1 core
#> ℹ [2026-01-22 03:00:59] Building results
#> ✔ [2026-01-22 03:00:59] Inferring network done
#> ℹ [2026-01-22 03:00:59] Network information:
#> ℹ                         Edges Regulators Targets
#> ℹ                       1   306         18      18
calculate_ji(
  network_table,
  example_ground_truth
)
#> $metrics
#>   Metric Value
#> 1     JI 0.514
#> 
```
