# Plot expression data

Plot expression data

## Usage

``` r
plot_scatter(
  data,
  smoothing_method = "lm",
  group_colors = RColorBrewer::brewer.pal(9, "Set1"),
  title_color = "black",
  title = NULL,
  col_title = NULL,
  row_title = NULL,
  legend_title = NULL,
  legend_position = "bottom",
  margins = "both",
  marginal_type = NULL,
  margins_size = 10,
  compute_correlation = TRUE,
  compute_correlation_method = "pearson",
  keep_aspect_ratio = TRUE,
  facet = FALSE,
  se = FALSE,
  pointdensity = TRUE
)
```

## Arguments

- data:

  Input data.

- smoothing_method:

  Method for smoothing curve, \`lm\` or \`loess\`.

- group_colors, title_color:

  Group and title colors.

- title, col_title, row_title, legend_title:

  Plot titles.

- legend_position:

  The position of legend.

- margins, marginal_type, margins_size:

  Marginal plot controls.

- compute_correlation, compute_correlation_method:

  Correlation controls.

- keep_aspect_ratio:

  Whether to use a 1:1 aspect ratio.

- facet:

  Whether to facet by group.

- se:

  Whether to show smoothing uncertainty.

- pointdensity:

  Whether to show point density.

## Value

A ggplot object.
