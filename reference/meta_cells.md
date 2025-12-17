# Detection of metacells from single-cell gene expression matrix

This function detects metacells from a single-cell gene expression
matrix using dimensionality reduction and clustering techniques.

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

  A gene expression matrix where rows represent genes and columns
  represent cells.

- genes_use:

  A character vector specifying genes to be used for PCA dimensionality
  reduction. Default is `NULL`.

- genes_exclude:

  A character vector specifying genes to be excluded from PCA
  computation. Default is `NULL`.

- var_genes_num:

  Number of most variable genes to select when `genes_use` is not
  provided. Default is `min(1000, nrow(matrix))`.

- gamma:

  Default is `10`. Coarse-graining parameter defining the target ratio
  of input cells to output metacells (e.g., gamma=10 yields
  approximately n/10 metacells for n input cells).

- knn_k:

  Default is `5`. Number of nearest neighbors for constructing the
  cell-cell similarity network.

- do_scale:

  Whether to standardize (center and scale) gene expression values
  before PCA. Default is `TRUE`.

- pc_num:

  Default is `25`. Number of principal components to retain for
  downstream analysis.

- fast_pca:

  Default is `TRUE`. Whether to use the faster
  [irlba](https://rdrr.io/pkg/irlba/man/irlba.html) algorithm instead of
  standard PCA. Recommended for large datasets.

- do_approx:

  Default is `FALSE`. Whether to use approximate nearest neighbor search
  for datasets with \>50000 cells to improve computational efficiency.

- approx_num:

  Default is `20000`. Number of cells to randomly sample for approximate
  nearest neighbor computation when `do_approx = TRUE`.

- directed:

  Default is `FALSE`. Whether to construct a directed or undirected
  nearest neighbor graph.

- use_nn2:

  Default is `TRUE`. Whether to use the faster RANN::nn2 algorithm for
  nearest neighbor search (only applicable with Euclidean distance).

- seed:

  Default is `1`. Random seed for reproducibility when subsampling cells
  in approximate mode.

- cluster_method:

  Default is `walktrap`. Algorithm for community detection in the cell
  similarity network. Options: `walktrap` (recommended) or `louvain`
  (gamma parameter ignored).

- block_size:

  Default is `10000`. Number of cells to process in each batch when
  mapping cells to metacells in approximate mode. Adjust based on
  available memory.

- weights:

  Default is `NULL`. Numeric vector of cell-specific weights for
  weighted averaging when computing metacell expression profiles. Length
  must match number of cells.

- do_median_norm:

  Default is `FALSE`. Whether to perform median-based normalization of
  the final metacell expression matrix.

- ...:

  Additional arguments passed to internal functions.

## Value

A matrix where rows represent metacells and columns represent genes.

## References

https://github.com/GfellerLab/SuperCell
https://github.com/kuijjerlab/SCORPION

## Examples

``` r
data(example_matrix)
meta_cells_matrix <- meta_cells(
  example_matrix
)
#> ! [2025-12-17 14:37:10] Number of PCs of PCA result is less than the desired number, using all PCs.
dim(meta_cells_matrix)
#> [1] 500  18
meta_cells_matrix[1:6, 1:6]
#> 6 x 6 Matrix of class "dgeMatrix"
#>            g1        g10        g11        g12        g13        g14
#> [1,] 2.180972 0.01930300 0.01449894 0.01160265 0.02435281 0.01841153
#> [2,] 2.372229 2.10205419 1.80734123 0.15135284 0.01745755 0.02049155
#> [3,] 2.300118 2.08998564 2.09137782 2.14921081 2.04650528 1.98375086
#> [4,] 2.042536 2.10596609 2.17743746 2.08665213 1.95181949 2.08697567
#> [5,] 2.666293 0.01307076 0.02208754 0.02237016 0.03820944 0.01280634
#> [6,] 1.965666 0.02490150 0.02487194 0.02585257 0.01419496 0.02565518
```
