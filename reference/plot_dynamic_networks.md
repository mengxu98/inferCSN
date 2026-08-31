# Plot dynamic networks

Plot dynamic networks

## Usage

``` r
plot_dynamic_networks(
  network_table,
  celltypes_order,
  ntop = 10,
  title = NULL,
  theme_type = "theme_void",
  plot_type = "ggplot",
  layout = "fruchtermanreingold",
  nrow = 2,
  figure_save = FALSE,
  figure_name = NULL,
  figure_width = 6,
  figure_height = 6,
  seed = 1
)
```

## Arguments

- network_table:

  Network edge table.

- celltypes_order:

  Cell-type order.

- ntop:

  Number of top genes to plot.

- title:

  Figure title.

- theme_type:

  Plot theme.

- plot_type:

  Output type.

- layout:

  Network layout.

- nrow:

  Number of rows.

- figure_save:

  Whether to save the figure.

- figure_name:

  Output filename.

- figure_width, figure_height:

  Figure dimensions.

- seed:

  Random seed.

## Value

A dynamic network plot.
