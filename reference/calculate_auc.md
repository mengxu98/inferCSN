# Calculate AUC Metrics

Calculates AUROC and AUPRC metrics with optional visualization

## Usage

``` r
calculate_auc(
  network_table,
  ground_truth,
  return_plot = FALSE,
  line_color = "#1563cc",
  line_width = 1,
  tag_levels = "A"
)
```

## Arguments

- network_table:

  A data frame of predicted network structure

- ground_truth:

  A data frame of ground truth network

- return_plot:

  Logical value indicating whether to generate plots

- line_color:

  Color for plot lines

- line_width:

  Width for plot lines

- tag_levels:

  Tag levels for plot annotations

## Value

A list containing metrics and optional plots

## Examples

``` r
data(example_matrix)
data("example_ground_truth")
network_table <- inferCSN(example_matrix)
#> ℹ [2026-01-09 07:09:25] Inferring network for <dense matrix>...
#> ◌ [2026-01-09 07:09:25] Checking parameters...
#> ℹ [2026-01-09 07:09:25] Using L0 sparse regression model
#> ℹ [2026-01-09 07:09:25] Using 1 core
#> ℹ [2026-01-09 07:09:25] Building results
#> ✔ [2026-01-09 07:09:25] Inferring network done
#> ℹ [2026-01-09 07:09:25] Network information:
#> ℹ                         Edges Regulators Targets
#> ℹ                       1   306         18      18
calculate_auc(
  network_table,
  example_ground_truth,
  return_plot = TRUE
)
#> $metrics
#>   Metric Value
#> 1  AUROC 0.952
#> 2  AUPRC 0.437
#> 
#> $plot

#> 
```
