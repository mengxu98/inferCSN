# Network Edge Comparison Visualization

Generates visualizations comparing edges of two networks.

## Usage

``` r
plot_edges_comparison(
  network_table,
  ground_truth,
  color_pattern = list(predicted = "gray", ground_truth = "#bb141a", overlap = "#1966ad",
    total = "#6C757D")
)
```

## Arguments

- network_table:

  A data frame of predicted network structure.

- ground_truth:

  A data frame of ground truth network.

- color_pattern:

  A list of colors for different categories.

## Value

A patchwork plot object containing network edge comparison and
distribution plots

## Examples

``` r
data(example_matrix)
data("example_ground_truth")
network_table <- inferCSN(example_matrix)
#> ℹ [2025-12-17 14:37:17] Inferring network for <dense matrix>...
#> ◌ [2025-12-17 14:37:17] Checking parameters...
#> ℹ [2025-12-17 14:37:17] Using "L0" sparse regression model
#> ℹ [2025-12-17 14:37:17] Using 1 core
#> ⠙ [2025-12-17 14:37:17] Running [1/18] Processing: g1  ETA:  0s
#> ✔ [2025-12-17 14:37:17] Completed 18 tasks in 183ms
#> 
#> ℹ [2025-12-17 14:37:17] Building results
#> ✔ [2025-12-17 14:37:18] Building network done
plot_edges_comparison(
  network_table,
  example_ground_truth
)
```
