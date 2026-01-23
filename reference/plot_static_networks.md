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
#> ℹ [2026-01-23 02:16:15] Inferring network for <matrix/array>...
#> ◌ [2026-01-23 02:16:15] Checking parameters...
#> ℹ [2026-01-23 02:16:15] Using L0 sparse regression model
#> ℹ [2026-01-23 02:16:15] Using 1 core
#> ⠙ [2026-01-23 02:16:15] Running for g1 [1/18] ■■■                              …
#> ✔ [2026-01-23 02:16:15] Completed 18 tasks in 187ms
#> 
#> ℹ [2026-01-23 02:16:15] Building results
#> ✔ [2026-01-23 02:16:16] Inferring network done
#> ℹ [2026-01-23 02:16:16] Network information:
#> ℹ                         Edges Regulators Targets
#> ℹ                       1   306         18      18
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
