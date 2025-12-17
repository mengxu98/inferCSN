# inferring cell-type specific gene regulatory network

inferring cell-type specific gene regulatory network

## Usage

``` r
inferCSN(
  object,
  penalty = "L0",
  cross_validation = FALSE,
  seed = 1,
  n_folds = 5,
  subsampling_method = c("sample", "meta_cells", "pseudobulk"),
  subsampling_ratio = 1,
  r_squared_threshold = 0,
  regulators = NULL,
  targets = NULL,
  cores = 1,
  verbose = TRUE,
  ...
)

# S4 method for class 'matrix'
inferCSN(
  object,
  penalty = "L0",
  cross_validation = FALSE,
  seed = 1,
  n_folds = 5,
  subsampling_method = c("sample", "meta_cells", "pseudobulk"),
  subsampling_ratio = 1,
  r_squared_threshold = 0,
  regulators = NULL,
  targets = NULL,
  cores = 1,
  verbose = TRUE,
  ...
)

# S4 method for class 'sparseMatrix'
inferCSN(
  object,
  penalty = "L0",
  cross_validation = FALSE,
  seed = 1,
  n_folds = 5,
  subsampling_method = c("sample", "meta_cells", "pseudobulk"),
  subsampling_ratio = 1,
  r_squared_threshold = 0,
  regulators = NULL,
  targets = NULL,
  cores = 1,
  verbose = TRUE,
  ...
)
```

## Arguments

- object:

  The input data for `inferCSN`.

- penalty:

  The type of regularization, default is `"L0"`. This can take either
  one of the following choices: `"L0"`, `"L0L1"`, and `"L0L2"`. For
  high-dimensional and sparse data, `"L0L2"` is more effective.

- cross_validation:

  Whether to use cross-validation. Default is `FALSE`.

- seed:

  The random seed for cross-validation. Default is `1`.

- n_folds:

  The number of folds for cross-validation. Default is `5`.

- subsampling_method:

  The method to use for subsampling. Options are `"sample"`,
  `"pseudobulk"` or `"meta_cells"`.

- subsampling_ratio:

  The percent of all samples used for
  [fit_srm](https://mengxu98.github.io/inferCSN/reference/fit_srm.md).
  Default is `1`.

- r_squared_threshold:

  Threshold of \\R^2\\ coefficient. Default is `0`.

- regulators:

  The regulator genes for which to infer the regulatory network.

- targets:

  The target genes for which to infer the regulatory network.

- cores:

  The number of cores to use for parallelization with
  [foreach::foreach](https://rdrr.io/pkg/foreach/man/foreach.html).
  Default is `1`.

- verbose:

  Whether to print progress messages. Default is `TRUE`.

- ...:

  Parameters for other methods.

## Value

A data table of regulator-target regulatory relationships. The data
table has the three columns: regulator, target, and weight.

## Examples

``` r
data(example_matrix)
network_table_1 <- inferCSN(
  example_matrix
)
#> ℹ [2025-12-17 14:37:09] Inferring network for <dense matrix>...
#> ◌ [2025-12-17 14:37:09] Checking parameters...
#> ℹ [2025-12-17 14:37:09] Using "L0" sparse regression model
#> ℹ [2025-12-17 14:37:09] Using 1 core
#> ⠙ [2025-12-17 14:37:09] Running [1/18] Processing: g1  ETA:  0s
#> ✔ [2025-12-17 14:37:09] Completed 18 tasks in 169ms
#> 
#> ℹ [2025-12-17 14:37:09] Building results
#> ✔ [2025-12-17 14:37:09] Building network done

network_table_2 <- inferCSN(
  example_matrix,
  cores = 2
)
#> ℹ [2025-12-17 14:37:09] Inferring network for <dense matrix>...
#> ◌ [2025-12-17 14:37:09] Checking parameters...
#> ℹ [2025-12-17 14:37:09] Using "L0" sparse regression model
#> ℹ [2025-12-17 14:37:09] Using 2 cores
#> ℹ [2025-12-17 14:37:09] Building results
#> ✔ [2025-12-17 14:37:09] Building network done

head(network_table_1)
#>   regulator target     weight
#> 1       g18     g1 -0.9223177
#> 2       g17    g18  0.8770468
#> 3        g4     g3  0.8103230
#> 4       g16    g15  0.7659245
#> 5       g17    g16  0.7558764
#> 6       g12    g11  0.7444053

identical(
  network_table_1,
  network_table_2
)
#> [1] TRUE

inferCSN(
  example_matrix,
  regulators = c("g1", "g2"),
  targets = c("g3", "g4")
)
#> ℹ [2025-12-17 14:37:09] Inferring network for <dense matrix>...
#> ◌ [2025-12-17 14:37:09] Checking parameters...
#> ℹ [2025-12-17 14:37:09] Using "L0" sparse regression model
#> ℹ [2025-12-17 14:37:09] Using 2 regulators
#> ℹ [2025-12-17 14:37:09] Using 2 targets
#> ℹ [2025-12-17 14:37:09] Using 1 core
#> ℹ [2025-12-17 14:37:09] Building results
#> ✔ [2025-12-17 14:37:09] Building network done
#>   regulator target     weight
#> 1        g2     g3  0.9848781
#> 2        g2     g4  0.9230387
#> 3        g1     g4 -0.3847071
#> 4        g1     g3 -0.1732490
inferCSN(
  example_matrix,
  regulators = c("g1", "g2"),
  targets = c("g3", "g0")
)
#> ℹ [2025-12-17 14:37:09] Inferring network for <dense matrix>...
#> ◌ [2025-12-17 14:37:09] Checking parameters...
#> ℹ [2025-12-17 14:37:09] Using "L0" sparse regression model
#> ℹ [2025-12-17 14:37:09] Using 2 regulators
#> ! [2025-12-17 14:37:09] 1 out of 2 targets are in the input matrix
#> ℹ [2025-12-17 14:37:09] Using 1 core
#> ℹ [2025-12-17 14:37:09] Building results
#> ✔ [2025-12-17 14:37:09] Building network done
#>   regulator target     weight
#> 1        g2     g3  0.9848781
#> 2        g1     g3 -0.1732490

if (FALSE) { # \dontrun{
data("example_ground_truth")
network_table_07 <- inferCSN(
  example_matrix,
  r_squared_threshold = 0.7
)
calculate_metrics(
  network_table_1,
  example_ground_truth,
  return_plot = TRUE
)
calculate_metrics(
  network_table_07,
  example_ground_truth,
  return_plot = TRUE
)
} # }
if (FALSE) { # \dontrun{
data(example_matrix)
network_table <- inferCSN(example_matrix)
head(network_table)

network_table_sparse_1 <- inferCSN(
  as(example_matrix, "sparseMatrix")
)
head(network_table_sparse_1)

network_table_sparse_2 <- inferCSN(
  as(example_matrix, "sparseMatrix"),
  cores = 2
)
identical(
  network_table,
  network_table_sparse_1
)

identical(
  network_table_sparse_1,
  network_table_sparse_2
)

plot_scatter(
  data.frame(
    network_table$weight,
    network_table_sparse_1$weight
  ),
  legend_position = "none"
)

plot_weight_distribution(
  network_table
) + plot_weight_distribution(
  network_table_sparse_1
)
} # }
```
