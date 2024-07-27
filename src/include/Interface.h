#ifndef RINTERFACE_H
#define RINTERFACE_H

#include "FitResult.h"
#include "Grid.h"
#include "GridParams.h"
#include "RcppArmadillo.h"
#include <string>
#include <vector>

inline void to_arma_error() {
  Rcpp::stop("model.fit only supports sparse matricies (dgCMatrix), 2D arrays "
             "(Dense Matricies)");
}

#endif // RINTERFACE_H
