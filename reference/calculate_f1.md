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
#> ℹ [2025-12-20 13:39:28] Inferring network for <dense matrix>...
#> ◌ [2025-12-20 13:39:28] Checking parameters...
#> ℹ [2025-12-20 13:39:28] Using "L0" sparse regression model
#> ℹ [2025-12-20 13:39:28] Using 1 core
#> ℹ [2025-12-20 13:39:28] Building results
#> ✔ [2025-12-20 13:39:28] Building network done
calculate_f1(
  network_table,
  example_ground_truth
)
#> $metrics
#>   Metric Value
#> 1     F1 0.692
#> 
```
