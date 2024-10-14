#include <Rcpp.h>
#include <algorithm>

using namespace Rcpp;

// [[Rcpp::export]]
DataFrame matrixToTable(NumericMatrix network_matrix) {
  int nrow = network_matrix.nrow();
  int ncol = network_matrix.ncol();

  CharacterVector row_names = rownames(network_matrix);
  CharacterVector col_names = colnames(network_matrix);

  CharacterVector regulators;
  CharacterVector targets;
  NumericVector weights;

  for (int i = 0; i < nrow; ++i) {
    for (int j = 0; j < ncol; ++j) {
      double weight = network_matrix(i, j);
      if (weight != 0) {
        regulators.push_back(row_names[i]);
        targets.push_back(col_names[j]);
        weights.push_back(weight);
      }
    }
  }

  return DataFrame::create(
    Named("regulator") = regulators,
    Named("target") = targets,
    Named("weight") = weights
  );
}
