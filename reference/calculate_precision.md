# Calculate Precision Metric

Calculates precision metric

## Usage

``` r
calculate_precision(network_table, ground_truth)
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
#> ℹ [2026-01-22 03:01:00] Inferring network for <dense matrix>...
#> ◌ [2026-01-22 03:01:00] Checking parameters...
#> ℹ [2026-01-22 03:01:00] Using L0 sparse regression model
#> ℹ [2026-01-22 03:01:00] Using 1 core
#> ℹ [2026-01-22 03:01:00] Building results
#> ✔ [2026-01-22 03:01:00] Inferring network done
#> ℹ [2026-01-22 03:01:00] Network information:
#> ℹ                         Edges Regulators Targets
#> ℹ                       1   306         18      18
calculate_precision(
  network_table,
  example_ground_truth
)
#> $metrics
#>      Metric Value
#> 1 Precision 0.529
#> 
```
