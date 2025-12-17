# Calculate F1 Score

Calculates F1 score

## Usage

``` r
calculate_f1(network_table, ground_truth)
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
#> ℹ [2025-12-17 14:37:04] Inferring network for <dense matrix>...
#> ◌ [2025-12-17 14:37:04] Checking parameters...
#> ℹ [2025-12-17 14:37:04] Using "L0" sparse regression model
#> ℹ [2025-12-17 14:37:04] Using 1 core
#> ℹ [2025-12-17 14:37:04] Building results
#> ✔ [2025-12-17 14:37:04] Building network done
calculate_f1(
  network_table,
  example_ground_truth
)
#> $metrics
#>   Metric Value
#> 1     F1 0.692
#> 
```
