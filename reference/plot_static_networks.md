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

  Network edge table.

- regulators, targets:

  Nodes to include.

- legend_position:

  The position of legend.

## Value

A ggplot2 object

## Examples

``` r
data(example_matrix)
network_table <- inferCSN(example_matrix)
#> ℹ [2026-08-31 02:40:02] Inferring network for <matrix/array>...
#> ◌ [2026-08-31 02:40:02] Checking parameters...
#> ✔ [2026-08-31 02:40:02] Inferring network done
#> ℹ [2026-08-31 02:40:02] Network information:
#> ℹ                         Edges Regulators Targets
#> ℹ                       1    12          6       6
example_edge <- network_table[1, ]
plot_static_networks(
  network_table,
  regulators = example_edge$regulator
)

plot_static_networks(
  network_table,
  targets = example_edge$target
)

plot_static_networks(
  network_table,
  regulators = example_edge$regulator,
  targets = example_edge$target
)
```
