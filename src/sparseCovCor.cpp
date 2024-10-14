#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// Calculate column means of a sparse matrix
vec sparse_col_means(const arma::sp_mat &X) {
  return vec(sum(X, 0).t()) / static_cast<double>(X.n_rows);
}

// Calculate covariance and correlation matrices for sparse matrices
// [[Rcpp::export]]
List sparseCovCor(const arma::sp_mat &x,
                  const Nullable<arma::sp_mat> &y_nullable = R_NilValue) {
  int n = x.n_rows;
  vec mu_x = sparse_col_means(x); // Calculate column means of x

  mat covmat;
  mat cormat;

  if (y_nullable.isNull()) {
    // Single matrix case, calculate covariance matrix of x
    covmat = (mat(x.t() * x) - n * (mu_x * mu_x.t())) / (n - 1);
    // Calculate standard deviation vector
    vec sdvec = sqrt(diagvec(covmat));

    // TODO: check if this is necessary
    // Avoid division by zero
    // Now it will be warning when compiling
    // sdvec = sdvec + (sdvec == 0) * 1e-8;

    // Calculate correlation matrix
    cormat = covmat / (sdvec * sdvec.t());
  } else {
    // Two-matrix case, calculate covariance matrix of x and y
    const arma::sp_mat &y = as<arma::sp_mat>(y_nullable);
    if (x.n_rows != y.n_rows) {
      stop("x and y should have the same number of rows");
    }

    vec mu_y = sparse_col_means(y); // Calculate column means of y

    // Calculate covariance matrix
    covmat = (mat(x.t() * y) - n * (mu_x * mu_y.t())) / (n - 1);

    // Calculate standard deviation vectors for x and y
    vec sdvec_x = sqrt((sum(square(x), 0).t() - n * square(mu_x)) / (n - 1));
    vec sdvec_y = sqrt((sum(square(y), 0).t() - n * square(mu_y)) / (n - 1));

    // TODO: check if this is necessary
    // Avoid division by zero
    // Now it will be warning when compiling
    // sdvec_x = sdvec_x + (sdvec_x == 0) * 1e-8;
    // sdvec_y = sdvec_y + (sdvec_y == 0) * 1e-8;

    // Calculate correlation matrix
    cormat = covmat / (sdvec_x * sdvec_y.t());
  }

  return List::create(Named("cov") = covmat, Named("cor") = cormat);
}

/*
# raw R code version of sparseCovCor

#' @title Fast correlation and covariance calcualtion for sparse matrices
#'
#' @inheritParams sparse_cor
sparseCovCor <- function(x, y = NULL) {
  if (!methods::is(x, "sparseMatrix")) {
    log_message(
      "x should be a sparse matrix",
      message_type = "error"
    )
  }
  n <- nrow(x)
  mu_x <- colMeans(x)
  if (is.null(y)) {
    covmat <- (as.matrix(crossprod(x)) - n * tcrossprod(mu_x)) / (n - 1)
    sdvec <- sqrt(diag(covmat))
    cormat <- covmat / tcrossprod(sdvec)
  } else {
    if (!methods::is(y, "sparseMatrix")) {
      log_message(
        "y should be a sparse matrix",
        message_type = "error"
      )
    }
    if (nrow(x) != nrow(y)) {
      log_message(
        "x and y should have the same number of rows",
        message_type = "error"
      )
    }

    mu_y <- colMeans(y)
    covmat <- (as.matrix(crossprod(x, y)) - n * tcrossprod(mu_x, mu_y)) / (n -
1) sdvecX <- sqrt((colSums(x^2) - n * mu_x^2) / (n - 1)) sdvecY <-
sqrt((colSums(y^2) - n * mu_y^2) / (n - 1)) cormat <- covmat /
tcrossprod(sdvecX, sdvecY)
  }
  return(
    list(
      cov = covmat,
      cor = cormat
    )
  )
}
*/
