% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/meta_cells.R
\name{meta_cells}
\alias{meta_cells}
\title{Detection of metacells with the SuperCell approach}
\usage{
meta_cells(
  matrix,
  genes_use = NULL,
  genes_exclude = NULL,
  var_genes_num = min(1000, nrow(matrix)),
  gamma = 10,
  knn_k = 5,
  do_scale = TRUE,
  pc_num = 25,
  fast_pca = TRUE,
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
}
\arguments{
\item{matrix}{log-normalized gene expression matrix with rows to be genes and cols to be cells}

\item{genes_use}{a vector of genes used to compute PCA}

\item{genes_exclude}{a vector of genes to be excluded when computing PCA}

\item{var_genes_num}{if \code{"genes.use"} is not provided, \code{"var_genes_num"} genes with the largest variation are used}

\item{gamma}{graining level of data (proportion of number of single cells in the initial dataset to the number of metacells in the final dataset)}

\item{knn_k}{parameter to compute single-cell kNN network}

\item{do_scale}{whether to scale gene expression matrix when computing PCA}

\item{pc_num}{number of principal components to use for construction of single-cell kNN network}

\item{fast_pca}{use \link[irlba]{irlba} as a faster version of prcomp}

\item{do_approx}{compute approximate kNN in case of a large dataset (>50'000)}

\item{approx_num}{number of cells to subsample for an approximate approach}

\item{directed}{directed}

\item{use_nn2}{use_nn2}

\item{seed}{seed to use to subsample cells for an approximate approach}

\item{cluster_method}{clustering method to identify metacells (available methods "walktrap" (default) and "louvain" (not recommended, gamma is ignored)).}

\item{block_size}{number of cells to map to the nearest metacell at the time (for approx coarse-graining)}

\item{weights}{vector of a cell weight (NULL by default), used for computing average gene expression withing cluster of metaells}

\item{do_median_norm}{whether to normalize by median value (FALSE by default)}

\item{...}{Parameters for other methods.}
}
\value{
a metacell matrix
}
\description{
This function detects metacells (former super-cells) from single-cell gene expression matrix
}
\examples{
data("example_matrix")
meta_cells_matrix <- meta_cells(
  example_matrix,
  pc_num = 10,
  fast_pca = FALSE
)
dim(meta_cells_matrix)
meta_cells_matrix[1:6,1:6]
}
