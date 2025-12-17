# Calculate Recall Metric

Calculates recall metric

## Usage

``` r
calculate_recall(network_table, ground_truth)
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
#> ℹ [2025-12-17 14:37:06] Inferring network for <dense matrix>...
#> ◌ [2025-12-17 14:37:06] Checking parameters...
#> ℹ [2025-12-17 14:37:06] Using "L0" sparse regression model
#> ℹ [2025-12-17 14:37:06] Using 1 core
#> ℹ [2025-12-17 14:37:06] Building results
#> ✔ [2025-12-17 14:37:07] Building network done
calculate_recall(
  network_table,
  example_ground_truth
)
#> $metrics
#>   Metric Value
#> 1 Recall     1
#> 
```
