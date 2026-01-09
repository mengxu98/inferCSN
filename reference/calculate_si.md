# Calculate Set Intersection

Calculates Set Intersection (SI) metric

## Usage

``` r
calculate_si(network_table, ground_truth)
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
#> ℹ [2026-01-09 07:09:30] Inferring network for <dense matrix>...
#> ◌ [2026-01-09 07:09:30] Checking parameters...
#> ℹ [2026-01-09 07:09:30] Using L0 sparse regression model
#> ℹ [2026-01-09 07:09:30] Using 1 core
#> ℹ [2026-01-09 07:09:30] Building results
#> ✔ [2026-01-09 07:09:30] Inferring network done
#> ℹ [2026-01-09 07:09:30] Network information:
#> ℹ                         Edges Regulators Targets
#> ℹ                       1   306         18      18
calculate_si(
  network_table,
  example_ground_truth
)
#> $metrics
#>   Metric Value
#> 1     SI    18
#> 
```
