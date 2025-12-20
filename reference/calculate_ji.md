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
#> ℹ [2025-12-20 13:39:29] Inferring network for <dense matrix>...
#> ◌ [2025-12-20 13:39:29] Checking parameters...
#> ℹ [2025-12-20 13:39:29] Using "L0" sparse regression model
#> ℹ [2025-12-20 13:39:29] Using 1 core
#> ⠙ [2025-12-20 13:39:29] Running [1/18] Processing: g1  ETA:  0s
#> ✔ [2025-12-20 13:39:29] Completed 18 tasks in 169ms
#> 
#> ℹ [2025-12-20 13:39:29] Building results
#> ✔ [2025-12-20 13:39:29] Building network done
calculate_ji(
  network_table,
  example_ground_truth
)
#> $metrics
#>   Metric Value
#> 1     JI 0.514
#> 
```
