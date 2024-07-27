#include <Rcpp.h>
#include <RcppParallel.h>
#include <algorithm>
#include <cmath>

using namespace Rcpp;
using namespace RcppParallel;

// [[Rcpp::export]]
NumericMatrix asMatrix(NumericVector rp, NumericVector cp, NumericVector z,
                       int nrows, int ncols) {
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
      std::string errMsg = "Index out of bounds at position (" +
                           std::to_string(rp[i]) + ", " +
                           std::to_string(cp[i]) + ").";
      stop(errMsg.c_str());
    }
    matrix(rp[i], cp[i]) = z[i];
  }

  return matrix;
}

// Parallel version reference: https://rcppcore.github.io/RcppParallel/
// [[Rcpp::depends(RcppParallel)]]
struct MatrixFiller : public Worker {
  const RVector<double> rp;
  const RVector<double> cp;
  const RVector<double> z;
  RMatrix<double> outputMatrix;

  MatrixFiller(const NumericVector rp, const NumericVector cp,
               const NumericVector z, NumericMatrix outputMatrix)
      : rp(rp), cp(cp), z(z), outputMatrix(outputMatrix) {}

  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end; ++i) {
      outputMatrix(rp[i], cp[i]) = z[i];
    }
  }
};

// [[Rcpp::export]]
NumericMatrix asMatrixParallel(NumericVector rp, NumericVector cp,
                               NumericVector z, int nrows, int ncols) {
  NumericMatrix outputMatrix(nrows, ncols);

  MatrixFiller matrixFiller(rp, cp, z, outputMatrix);

  int grainSize = std::max(1, static_cast<int>(z.size() / 2000));
  RcppParallel::parallelFor(0, z.size(), matrixFiller, grainSize);

  return outputMatrix;
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
 i <- sample(1:2000, 500) # Row indices (assuming we have 50 non-zero elements)
 j <- sample(1:2000, 500) # Column indices
 x <- rnorm(100) # Corresponding non-zero values
 dimnames <- list(
 paste0("a", rep(1:2000)),
 paste0("b", rep(1:2000))
 )

# Create the sparse matrix
 sparse_matrix <- Matrix::sparseMatrix(
 i = i,
 j = j,
 x = x,
 dims = c(2000, 2000),
 dimnames = dimnames
 )

 identical(
 as.matrix(sparse_matrix),
 as_matrix(sparse_matrix)
 )

 identical(
 as.matrix(sparse_matrix),
 as_matrix(sparse_matrix, parallel = TRUE)
 )
 if (!requireNamespace("bench", quietly = TRUE)) {
 pak::pak('bench')
 }
 bench::mark(
  as.matrix(sparse_matrix),
  as_matrix(sparse_matrix),
  as_matrix(sparse_matrix, parallel = TRUE)
 )

 if (!requireNamespace("rbenchmark", quietly = TRUE)) {
 pak::pak('rbenchmark')
 }
 rbenchmark::benchmark(
 as.matrix(sparse_matrix),
 as_matrix(sparse_matrix),
 as_matrix(sparse_matrix, parallel = TRUE)
 )

 */
