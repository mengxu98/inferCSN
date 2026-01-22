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
#> ℹ [2026-01-22 03:00:58] Inferring network for <dense matrix>...
#> ◌ [2026-01-22 03:00:58] Checking parameters...
#> ℹ [2026-01-22 03:00:58] Using L0 sparse regression model
#> ℹ [2026-01-22 03:00:58] Using 1 core
#> ⠙ [2026-01-22 03:00:58] Running for g1 [1/18] ■■■                              …
#> ✔ [2026-01-22 03:00:58] Completed 18 tasks in 212ms
#> 
#> ℹ [2026-01-22 03:00:58] Building results
#> ✔ [2026-01-22 03:00:58] Inferring network done
#> ℹ [2026-01-22 03:00:58] Network information:
#> ℹ                         Edges Regulators Targets
#> ℹ                       1   306         18      18
calculate_f1(
  network_table,
  example_ground_truth
)
#> $metrics
#>   Metric Value
#> 1     F1 0.692
#> 
```
