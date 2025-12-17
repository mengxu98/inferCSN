# Plot dynamic networks

Plot dynamic networks

## Usage

``` r
plot_static_networks(
  network_table,
  regulators = NULL,
  targets = NULL,
  legend_position = "right"
)
```

## Arguments

- network_table:

  The weight data table of network.

- regulators:

  Regulators list.

- targets:

  Targets list.

- legend_position:

  The position of legend.

## Value

A ggplot2 object

## Examples

``` r
data(example_matrix)
network_table <- inferCSN(example_matrix)
#> ℹ [2025-12-17 14:37:29] Inferring network for <dense matrix>...
#> ◌ [2025-12-17 14:37:29] Checking parameters...
#> ℹ [2025-12-17 14:37:29] Using "L0" sparse regression model
#> ℹ [2025-12-17 14:37:29] Using 1 core
#> ⠙ [2025-12-17 14:37:29] Running [1/18] Processing: g1  ETA:  0s
#> ✔ [2025-12-17 14:37:29] Completed 18 tasks in 195ms
#> 
#> ℹ [2025-12-17 14:37:29] Building results
#> ✔ [2025-12-17 14:37:29] Building network done
plot_static_networks(
  network_table,
  regulators = "g1"
)

plot_static_networks(
  network_table,
  targets = "g1"
)

plot_static_networks(
  network_table,
  regulators = "g1",
  targets = "g2"
)
```
