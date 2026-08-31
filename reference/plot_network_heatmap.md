# Plot network heatmaps

Plot network heatmaps

## Usage

``` r
plot_network_heatmap(
  network_table,
  regulators = NULL,
  targets = NULL,
  switch_matrix = TRUE,
  show_names = FALSE,
  show_names_position = c("outer", "all"),
  rect_color = NA,
  heatmap_size_lock = TRUE,
  heatmap_size = 5,
  heatmap_height = NULL,
  heatmap_width = NULL,
  heatmap_title = NULL,
  ncol = NULL,
  nrow = NULL,
  performance_metrics = "none",
  performance_ground_truth = NULL,
  truth_cell_border = TRUE,
  truth_cell_border_legend = TRUE,
  truth_match_border_color = "#2E7D32",
  truth_mismatch_border_color = "#F2C94C",
  truth_cell_border_lwd = 1.2,
  border_color = "gray",
  heatmap_palette = "viridis",
  heatmap_palcolor = NULL,
  row_anno_palette = "Set1",
  row_anno_palcolor = NULL,
  col_anno_palette = "Set2",
  col_anno_palcolor = NULL,
  row_anno = FALSE,
  column_anno = FALSE,
  anno_width = 1,
  anno_height = 1,
  row_anno_type = c("boxplot", "barplot", "histogram", "density", "lines", "points",
    "horizon"),
  column_anno_type = c("boxplot", "barplot", "histogram", "density", "lines", "points"),
  legend_name = NULL,
  row_title = "Regulators"
)
```

## Arguments

- network_table:

  Network edge table.

- regulators, targets:

  Nodes to include.

- switch_matrix:

  Convert edge tables to matrices.

- show_names, heatmap_size_lock:

  Logical display controls.

- show_names_position:

  Name placement for multiple heatmaps.

- heatmap_size, heatmap_height, heatmap_width:

  Heatmap dimensions.

- heatmap_title:

  Heatmap title.

- ncol, nrow:

  Layout dimensions for multiple heatmaps.

- performance_metrics:

  Metrics appended to titles.

- performance_ground_truth:

  Ground-truth network for metrics.

- truth_cell_border, truth_cell_border_legend:

  Ground-truth border controls.

- truth_match_border_color, truth_mismatch_border_color:

  Border colors.

- truth_cell_border_lwd:

  Border width.

- border_color, rect_color:

  Border and cell colors.

- heatmap_palcolor, heatmap_palette:

  Heatmap colors.

- row_anno_palette, row_anno_palcolor:

  Row annotation colors.

- col_anno_palette, col_anno_palcolor:

  Column annotation colors.

- row_anno, column_anno:

  Enable annotations.

- anno_width, anno_height:

  Annotation dimensions.

- row_anno_type, column_anno_type:

  Annotation plot types.

- legend_name, row_title:

  Legend and row titles.

## Value

A heatmap object.
