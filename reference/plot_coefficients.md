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
#> ℹ [2025-12-20 13:39:35] Inferring network for <dense matrix>...
#> ◌ [2025-12-20 13:39:35] Checking parameters...
#> ℹ [2025-12-20 13:39:35] Using "L0" sparse regression model
#> ℹ [2025-12-20 13:39:35] Using 3 targets
#> ℹ [2025-12-20 13:39:35] Using 1 core
#> ⠙ [2025-12-20 13:39:35] Running [2/3] Processing: g2  ETA:  0s
#> ✔ [2025-12-20 13:39:35] Completed 3 tasks in 47ms
#> 
#> ℹ [2025-12-20 13:39:35] Building results
#> ✔ [2025-12-20 13:39:35] Building network done
plot_coefficients(network_table, show_values = FALSE)

plot_coefficients(network_table, targets = "g1")
```
