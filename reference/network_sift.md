# Sifting network

Sifting network

## Usage

``` r
network_sift(
  network_table,
  matrix = NULL,
  meta_data = NULL,
  pseudotime_column = NULL,
  method = c("entropy", "max"),
  entropy_method = c("Shannon", "Renyi"),
  effective_entropy = FALSE,
  shuffles = 100,
  entropy_nboot = 300,
  lag_value = 1,
  entropy_p_value = 0.05,
  cores = 1,
  verbose = TRUE
)
```

## Arguments

- network_table:

  The weight data table of network.

- matrix:

  The expression matrix.

- meta_data:

  The meta data for cells or samples.

- pseudotime_column:

  The column of pseudotime.

- method:

  The method used for filter edges. Could be choose `"entropy"` or
  `"max"`.

- entropy_method:

  If setting `method` to `"entropy"`, could be choose `"Shannon"` or
  `"Renyi"` to compute entropy.

- effective_entropy:

  Default is `FALSE`. Whether to use effective entropy to filter
  weights.

- shuffles:

  Default is `100`. The number of shuffles used to calculate the
  effective transfer entropy.

- entropy_nboot:

  Default is `300`. The number of bootstrap replications for each
  direction of the estimated transfer entropy.

- lag_value:

  Default is `1`. Markov order of gene expression values, i.e. the
  number of lagged values affecting the current value of gene expression
  values.

- entropy_p_value:

  P value used to filter edges by entropy. Default is `0.05`.

- cores:

  The number of cores to use for parallelization with
  [foreach::foreach](https://rdrr.io/pkg/foreach/man/foreach.html).
  Default is `1`.

- verbose:

  Whether to print progress messages. Default is `TRUE`.

## Value

A data table of regulator-target regulatory relationships. The data
table has the three columns: regulator, target, and weight.

## Examples

``` r
if (FALSE) { # \dontrun{
data(example_matrix)
data(example_meta_data)
data(example_ground_truth)

network_table <- inferCSN(example_matrix)
network_table_sifted <- network_sift(network_table)
network_table_sifted_entropy <- network_sift(
  network_table,
  matrix = example_matrix,
  meta_data = example_meta_data,
  pseudotime_column = "pseudotime",
  lag_value = 2,
  shuffles = 0,
  entropy_nboot = 0
)

plot_network_heatmap(
  example_ground_truth[, 1:3],
  heatmap_title = "Ground truth",
  show_names = TRUE,
  rect_color = "gray70"
)
plot_network_heatmap(
  network_table,
  heatmap_title = "Raw",
  show_names = TRUE,
  rect_color = "gray70"
)
plot_network_heatmap(
  network_table_sifted,
  heatmap_title = "Filtered",
  show_names = TRUE,
  rect_color = "gray70"
)
plot_network_heatmap(
  network_table_sifted_entropy,
  heatmap_title = "Filtered by entropy",
  show_names = TRUE,
  rect_color = "gray70"
)

calculate_metrics(
  network_table,
  example_ground_truth,
  return_plot = TRUE
)
calculate_metrics(
  network_table_sifted,
  example_ground_truth,
  return_plot = TRUE
)
calculate_metrics(
  network_table_sifted_entropy,
  example_ground_truth,
  return_plot = TRUE
)
} # }
```
