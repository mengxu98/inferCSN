# Calculate Accuracy

Calculates accuracy metric

## Usage

``` r
calculate_accuracy(network_table, ground_truth)
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
#> ℹ [2025-12-20 13:39:26] Inferring network for <dense matrix>...
#> ◌ [2025-12-20 13:39:26] Checking parameters...
#> ℹ [2025-12-20 13:39:26] Using "L0" sparse regression model
#> ℹ [2025-12-20 13:39:26] Using 1 core
#> ⠙ [2025-12-20 13:39:26] Running [1/18] Processing: g1  ETA:  1s
#> ✔ [2025-12-20 13:39:26] Completed 18 tasks in 218ms
#> 
#> ℹ [2025-12-20 13:39:26] Building results
#> ✔ [2025-12-20 13:39:26] Building network done
calculate_accuracy(
  network_table,
  example_ground_truth
)
#> $metrics
#>   Metric Value
#> 1    ACC 0.948
#> 
```
