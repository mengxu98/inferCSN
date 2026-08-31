# Plot coefficients for multiple targets

Plot coefficients for multiple targets

## Usage

``` r
plot_coefficients(data, targets = NULL, nrow = NULL, ...)
```

## Arguments

- data:

  Input data.

- targets:

  Targets to plot. Default is \`NULL\`.

- nrow:

  Number of rows for the plot. Default is \`NULL\`.

- ...:

  Other arguments passed to \[plot_coefficient\].

## Value

A list of ggplot objects

## Examples

``` r
data(example_matrix)
network_table <- inferCSN(
  example_matrix,
  targets = c("g1", "g2", "g3")
)
#> ℹ [2026-08-31 02:39:59] Inferring network for <matrix/array>...
#> ◌ [2026-08-31 02:39:59] Checking parameters...
#> ✔ [2026-08-31 02:39:59] Inferring network done
#> ℹ [2026-08-31 02:39:59] Network information:
#> ℹ                         Edges Regulators Targets
#> ℹ                       1     6          5       3
plot_coefficients(network_table, show_values = FALSE)

plot_coefficients(network_table, targets = "g1")
```
