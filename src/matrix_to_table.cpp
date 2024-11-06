#include <Rcpp.h>
#include <algorithm>
#include "network_format.h"

using namespace Rcpp;

//' @title Switch matrix to network table
//'
//' @inheritParams table_to_matrix
//' @param network_matrix The matrix of network weight.
//'
//' @return Network table
//' @export
//'
//' @examples
//' data("example_matrix")
//' network_table <- inferCSN(example_matrix)
//' network_matrix <- table_to_matrix(network_table)
//' network_table_new <- matrix_to_table(network_matrix)
//' head(network_table)
//' head(network_table_new)
//' identical(
//'   network_table,
//'   network_table_new
//' )
// [[Rcpp::export]]
DataFrame matrix_to_table(NumericMatrix network_matrix,
                          Nullable<CharacterVector> regulators = R_NilValue,
                          Nullable<CharacterVector> targets = R_NilValue)
{
  int nrow = network_matrix.nrow();
  int ncol = network_matrix.ncol();

  CharacterVector row_names = rownames(network_matrix);
  CharacterVector col_names = colnames(network_matrix);

  if (row_names.length() == 0 || col_names.length() == 0)
  {
    stop("Input matrix must have both row and column names");
  }

  CharacterVector regulators_out;
  CharacterVector targets_out;
  NumericVector weights;

  for (int i = 0; i < nrow; ++i)
  {
    for (int j = 0; j < ncol; ++j)
    {
      double weight = network_matrix(i, j);
      if (weight != 0)
      {
        regulators_out.push_back(row_names[i]);
        targets_out.push_back(col_names[j]);
        weights.push_back(weight);
      }
    }
  }

  DataFrame intermediate = DataFrame::create(
      Rcpp::Named("regulator") = regulators_out,
      Rcpp::Named("target") = targets_out,
      Rcpp::Named("weight") = weights);

  return network_format(intermediate, regulators, targets, false);
}

/*
matrix_to_table <- function(
    network_matrix,
    regulators = NULL,
    targets = NULL) {
  filter_sort_matrix(
    network_matrix,
    regulators = regulators,
    targets = targets
  ) |>
    matrixToTable() |>
    network_format(abs_weight = FALSE)
}
*/