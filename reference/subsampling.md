# Subsampling function

This function subsamples a matrix using either random sampling or meta
cells method.

## Usage

``` r
subsampling(
  matrix,
  subsampling_method = c("sample", "meta_cells", "pseudobulk"),
  subsampling_ratio = 1,
  seed = 1,
  verbose = TRUE,
  ...
)
```

## Arguments

- matrix:

  The input matrix to be subsampled.

- subsampling_method:

  The method to use for subsampling. Options are `"sample"`,
  `"pseudobulk"` or `"meta_cells"`.

- subsampling_ratio:

  The percent of all samples used for
  [fit_srm](https://mengxu98.github.io/inferCSN/reference/fit_srm.md).
  Default is `1`.

- seed:

  The random seed for cross-validation. Default is `1`.

- verbose:

  Whether to print progress messages. Default is `TRUE`.

- ...:

  Parameters for other methods.

## Value

The subsampled matrix.

## Examples

``` r
data(example_matrix)
data("example_ground_truth")
subsample_matrix <- subsampling(
  example_matrix,
  subsampling_ratio = 0.5
)
#> ℹ [2025-12-17 14:37:32] Subsample matrix generated, dimensions: 2500 cells by 18 genes
subsample_matrix_2 <- subsampling(
  example_matrix,
  subsampling_method = "meta_cells",
  subsampling_ratio = 0.5,
  fast_pca = FALSE
)
#> ! [2025-12-17 14:37:32] Number of PCs of PCA result is less than the desired number, using all PCs.
#> ℹ [2025-12-17 14:37:33] Subsample matrix generated, dimensions: 2500 cells by 18 genes
subsample_matrix_3 <- subsampling(
  example_matrix,
  subsampling_method = "pseudobulk",
  subsampling_ratio = 0.5
)
#> ℹ [2025-12-17 14:37:34] Subsample matrix generated, dimensions: 2500 cells by 18 genes

calculate_metrics(
  inferCSN(example_matrix),
  example_ground_truth,
  return_plot = TRUE
)
#> ℹ [2025-12-17 14:37:34] Inferring network for <dense matrix>...
#> ◌ [2025-12-17 14:37:34] Checking parameters...
#> ℹ [2025-12-17 14:37:34] Using "L0" sparse regression model
#> ℹ [2025-12-17 14:37:34] Using 1 core
#> ⠙ [2025-12-17 14:37:34] Running [1/18] Processing: g1  ETA:  0s
#> ✔ [2025-12-17 14:37:34] Completed 18 tasks in 195ms
#> 
#> ℹ [2025-12-17 14:37:34] Building results
#> ✔ [2025-12-17 14:37:34] Building network done
#> $metrics
#>      Metric  Value
#> 1     AUROC  0.952
#> 2     AUPRC  0.437
#> 3 Precision  0.529
#> 4    Recall  1.000
#> 5        F1  0.692
#> 6       ACC  0.948
#> 7        JI  0.514
#> 8        SI 18.000
#> 
#> $plot

#> 
calculate_metrics(
  inferCSN(subsample_matrix),
  example_ground_truth,
  return_plot = TRUE
)
#> ℹ [2025-12-17 14:37:34] Inferring network for <dense matrix>...
#> ◌ [2025-12-17 14:37:34] Checking parameters...
#> ℹ [2025-12-17 14:37:34] Using "L0" sparse regression model
#> ℹ [2025-12-17 14:37:34] Using 1 core
#> ℹ [2025-12-17 14:37:34] Building results
#> ✔ [2025-12-17 14:37:34] Building network done
#> $metrics
#>      Metric  Value
#> 1     AUROC  0.955
#> 2     AUPRC  0.449
#> 3 Precision  0.529
#> 4    Recall  1.000
#> 5        F1  0.692
#> 6       ACC  0.948
#> 7        JI  0.514
#> 8        SI 18.000
#> 
#> $plot

#> 
calculate_metrics(
  inferCSN(subsample_matrix_2),
  example_ground_truth,
  return_plot = TRUE
)
#> ℹ [2025-12-17 14:37:34] Inferring network for <dense matrix>...
#> ◌ [2025-12-17 14:37:34] Checking parameters...
#> ℹ [2025-12-17 14:37:34] Using "L0" sparse regression model
#> ℹ [2025-12-17 14:37:34] Using 1 core
#> ℹ [2025-12-17 14:37:34] Building results
#> ✔ [2025-12-17 14:37:34] Building network done
#> $metrics
#>      Metric  Value
#> 1     AUROC  0.952
#> 2     AUPRC  0.439
#> 3 Precision  0.529
#> 4    Recall  1.000
#> 5        F1  0.692
#> 6       ACC  0.948
#> 7        JI  0.514
#> 8        SI 18.000
#> 
#> $plot

#> 
calculate_metrics(
  inferCSN(subsample_matrix_3),
  example_ground_truth,
  return_plot = TRUE
)
#> ℹ [2025-12-17 14:37:35] Inferring network for <dense matrix>...
#> ◌ [2025-12-17 14:37:35] Checking parameters...
#> ℹ [2025-12-17 14:37:35] Using "L0" sparse regression model
#> ℹ [2025-12-17 14:37:35] Using 1 core
#> ℹ [2025-12-17 14:37:35] Building results
#> ✔ [2025-12-17 14:37:35] Building network done
#> $metrics
#>      Metric  Value
#> 1     AUROC  0.955
#> 2     AUPRC  0.449
#> 3 Precision  0.529
#> 4    Recall  1.000
#> 5        F1  0.692
#> 6       ACC  0.948
#> 7        JI  0.514
#> 8        SI 18.000
#> 
#> $plot

#> 
```
