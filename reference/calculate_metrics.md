# Calculate network metrics

Calculate network metrics

## Usage

``` r
calculate_metrics(
  network_table,
  ground_truth,
  metric_type = c("all", "all_no_epr", "auc", "auroc", "auprc", "epr", "precision",
    "recall", "f1", "accuracy", "si", "ji"),
  return_plot = FALSE,
  tf_edges = FALSE,
  line_color = "#1563cc",
  line_width = 1
)
```

## Arguments

- network_table:

  Predicted edge table.

- ground_truth:

  Ground-truth edge table.

- metric_type:

  Metrics to compute.

- return_plot:

  Whether to return a plot.

- tf_edges:

  Whether to restrict regulators to transcription factors.

- line_color, line_width:

  Plot line style.

## Value

A list with metrics and an optional plot.
