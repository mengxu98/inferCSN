#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix asMatrix(NumericVector rp,
                       NumericVector cp,
                       NumericVector z,
                       int nrows,
                       int ncols) {
  int elementsCount = z.size();

  // Input validation
  if (elementsCount != rp.size() || elementsCount != cp.size()) {
    stop("The lengths of 'rp', 'cp', and 'z' must be equal.");
  }
  if (nrows <= 0 || ncols <= 0) {
    stop("Both 'nrows' and 'ncols' must be positive.");
  }

  NumericMatrix matrix(nrows, ncols);

  // Index validity check and matrix population
  for (int i = 0; i < elementsCount; i++) {
    if (rp[i] < 0 || rp[i] >= nrows || cp[i] < 0 || cp[i] >= ncols) {
      std::string errMsg = "Index out of bounds at position (" + std::to_string(rp[i]) + ", " + std::to_string(cp[i]) + ").";
      stop(errMsg.c_str());
    }
    matrix(rp[i], cp[i]) = z[i];
  }

  return matrix;
}


// Reference:
// https://github.com/zhanghao-njmu/SCP/blob/b9b0eb7a7bf2c2c4b2262e73e09d7ebd515c7da0/R/utils.R#L831
// https://github.com/zhanghao-njmu/SCP/blob/b9b0eb7a7bf2c2c4b2262e73e09d7ebd515c7da0/src/asMatrix.cpp#L5
// IntegerMatrix asMatrix(NumericVector rp,
//                        NumericVector cp,
//                        NumericVector x,
//                        int nrows,
//                        int ncols) {
//   int k = x.size();
//   IntegerMatrix mat(nrows, ncols);
//   for (int i = 0; i < k; i++) {
//     mat(rp[i], cp[i]) = x[i];
//   }
//   return mat;
// }

/*

 Rcpp::sourceCpp("functions/asMatrix.cpp")

 as_matrix <- function(x) {
 if (!inherits(matrix, "dgCMatrix")) {
 return(as.matrix(x))
 } else {
 row_pos <- x@i
 col_pos <- findInterval(seq_along(x@x) - 1, x@p[-1])
 out <- asMatrix(
 rp = row_pos,
 cp = col_pos,
 z = x@x,
 nrows = x@Dim[1],
 ncols = x@Dim[2]
 )
 attr(out, "dimnames") <- list(x@Dimnames[[1]], x@Dimnames[[2]])
 return(out)
 }
 }

 if (!requireNamespace("Matrix", quietly = TRUE)) {
 install.packages("Matrix")
 }

# Define row indices, column indices, and corresponding non-zero values
 i <- sample(1:200, 50) # Row indices (assuming we have 50 non-zero elements)
 j <- sample(1:200, 50) # Column indices
 x <- rnorm(50) # Corresponding non-zero values
 dimnames <- list(
 paste0("a", rep(1:200)),
 paste0("b", rep(1:200))
 )

# Create the sparse matrix
 sparse_matrix <- Matrix::sparseMatrix(
 i = i,
 j = j,
 x = x,
 dims = c(200, 200),
 dimnames = dimnames
 )

 identical(
 as.matrix(sparse_matrix),
 as_matrix(sparse_matrix)
 )

 */
