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
#> ℹ [2026-01-23 02:15:49] Inferring network for <matrix/array>...
#> ◌ [2026-01-23 02:15:49] Checking parameters...
#> ℹ [2026-01-23 02:15:49] Using L0 sparse regression model
#> ℹ [2026-01-23 02:15:49] Using 1 core
#> ⠙ [2026-01-23 02:15:49] Running for g1 [1/18] ■■■                              …
#> ⠹ [2026-01-23 02:15:49] Running for g2 [11/18] ■■■■■■■■■■■■■■■■■■■             …
#> ✔ [2026-01-23 02:15:49] Completed 18 tasks in 267ms
#> 
#> ℹ [2026-01-23 02:15:49] Building results
#> ✔ [2026-01-23 02:15:50] Inferring network done
#> ℹ [2026-01-23 02:15:50] Network information:
#> ℹ                         Edges Regulators Targets
#> ℹ                       1   306         18      18
calculate_accuracy(
  network_table,
  example_ground_truth
)
#> $metrics
#>   Metric Value
#> 1    ACC 0.948
#> 
```
