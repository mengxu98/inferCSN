# Plot contrast networks

Plot contrast networks

## Usage

``` r
plot_contrast_networks(
  network_table,
  degree_value = 0,
  weight_value = 0,
  legend_position = "bottom"
)
```

## Arguments

- network_table:

  The weight data table of network.

- degree_value:

  Degree value to filter nodes. Default is `0`.

- weight_value:

  Weight value to filter edges. Default is `0`.

- legend_position:

  The position of legend. Default is `"bottom"`.

## Value

A ggplot2 object.

## Examples

``` r
data(example_matrix)
network_table <- inferCSN(example_matrix)
#> ℹ [2025-12-20 13:39:36] Inferring network for <dense matrix>...
#> ◌ [2025-12-20 13:39:36] Checking parameters...
#> ℹ [2025-12-20 13:39:36] Using "L0" sparse regression model
#> ℹ [2025-12-20 13:39:36] Using 1 core
#> ℹ [2025-12-20 13:39:36] Building results
#> ✔ [2025-12-20 13:39:36] Building network done
plot_contrast_networks(network_table[1:50, ])
```
