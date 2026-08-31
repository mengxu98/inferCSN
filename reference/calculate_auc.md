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
  tag_levels = "A",
  tf_edges = FALSE
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

- tf_edges:

  Whether to restrict candidate edges to TF-to-gene. Default is
  \`FALSE\`.

## Value

A list containing metrics and optional plots
