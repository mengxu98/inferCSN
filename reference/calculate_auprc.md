# Calculate AUPRC Metric

Calculates AUPRC metric with optional visualization

## Usage

``` r
calculate_auprc(
  network_table,
  ground_truth,
  return_plot = FALSE,
  line_color = "#1563cc",
  line_width = 1,
  tf_edges = FALSE
)
```

## Arguments

- network_table:

  A data frame of predicted network structure

- ground_truth:

  A data frame of ground truth network

- return_plot:

  Logical value indicating whether to generate plot

- line_color:

  Color for plot lines

- line_width:

  Width for plot lines

- tf_edges:

  Whether to restrict candidate edges to TF-to-gene. Default is
  \`FALSE\`.

## Value

A list containing metric and optional plot
