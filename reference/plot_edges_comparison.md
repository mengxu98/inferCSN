# Plot Edges Comparison

Creates a scatter plot comparing predicted, ground truth, and
overlapping edges between two gene regulatory networks.

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

  A data frame of predicted network structure with \`regulator\` and
  \`target\` columns.

- ground_truth:

  A data frame of ground truth network with \`regulator\` and \`target\`
  columns.

- color_pattern:

  A named list of colors for predicted, ground_truth, overlap, and total
  edges.

## Value

A ggplot object visualizing edge overlap between the two networks.
