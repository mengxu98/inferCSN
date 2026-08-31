# Build metacells

Build metacells

## Usage

``` r
meta_cells(
  matrix,
  genes_use = NULL,
  genes_exclude = NULL,
  var_genes_num = min(1000, nrow(matrix)),
  gamma = 10,
  knn_k = 5,
  do_scale = TRUE,
  pc_num = 25,
  fast_pca = FALSE,
  do_approx = FALSE,
  approx_num = 20000,
  directed = FALSE,
  use_nn2 = TRUE,
  seed = 1,
  cluster_method = "walktrap",
  block_size = 10000,
  weights = NULL,
  do_median_norm = FALSE,
  ...
)
```

## Arguments

- matrix:

  Cell-by-gene expression matrix.

- genes_use, genes_exclude:

  Genes to include or exclude from PCA.

- var_genes_num:

  Number of variable genes.

- gamma:

  Target cells-per-metacell ratio.

- knn_k:

  Number of nearest neighbors.

- do_scale:

  Whether to scale expression values.

- pc_num:

  Number of principal components.

- fast_pca:

  Whether to use truncated PCA.

- do_approx:

  Whether to use approximate neighbor search.

- approx_num:

  Number of cells used in approximate mode.

- directed:

  Whether to build a directed graph.

- use_nn2:

  Whether to use nearest-neighbor search.

- seed:

  Random seed.

- cluster_method:

  Community-detection method.

- block_size:

  Mapping batch size.

- weights:

  Optional cell weights.

- do_median_norm:

  Whether to median-normalize the result.

- ...:

  Additional arguments.

## Value

A metacell-by-gene matrix.
