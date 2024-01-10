#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix table_to_matrix(DataFrame weight_table) {
  CharacterVector regulator = weight_table[0];
  CharacterVector target = weight_table[1];
  NumericVector weight = weight_table[2];

  CharacterVector genes = union_(regulator, target);

  int num_genes = genes.size();

  NumericMatrix weight_matrix(num_genes, num_genes);
  colnames(weight_matrix) = genes;
  rownames(weight_matrix) = genes;

  for (int i = 0; i < weight_table.nrows(); ++i) {
    int regulators_index = std::distance(genes.begin(), std::find(genes.begin(), genes.end(), regulator[i]));
    int targets_index = std::distance(genes.begin(), std::find(genes.begin(), genes.end(), target[i]));
    weight_matrix(regulators_index, targets_index) = weight[i];
  }

  return weight_matrix;
}

// #' R code
// #' @title table_to_matrix_r
// #'  Switch weight data table to weight matrix
// #' @param weight_table The weight data table of network
// #'
// #' @return Weight matrix
// #' @export
// #'
// table_to_matrix_r <- function(weight_table) {
//   colnames(weight_table) <- c("regulator", "target", "weight")
//   genes <- gtools::mixedsort(unique(c(weight_table$regulator, weight_table$target)))
//
//   weight_matrix <- matrix(0,
//                           nrow = length(genes),
//                           ncol = length(genes),
//                           dimnames = list(genes, genes))
//
//   for (i in 1:nrow(weight_table)) {
//     weight_matrix[weight_table$regulator[i], weight_table$target[i]] <- weight_table$weight[i]
//   }
//
//   return(weight_matrix)
// }
//
// table_to_matrix_r2 <- function(weight_table) {
//   weight_matrix <- reshape2::acast(weight_table, regulator ~ target)[genes, genes]
//   weight_matrix[is.na(weight_matrix)] <- 0
//   return(weight_matrix)
// }
//
// #' @examples
//   Rcpp::sourceCpp("src/table_to_matrix.cpp")
//   library(inferCSN)
//   data("example_matrix")
//   weight_table <- inferCSN(example_matrix, verbose = TRUE)
//   weight_matrix <- table_to_matrix(weight_table)
//   genes <- gtools::mixedsort(unique(c(weight_table$regulator, weight_table$target)))
//   weight_matrix <- weight_matrix[genes, genes]
//   weight_matrix_r <- table_to_matrix_r(weight_table)
//   weightMatrixR2 <- table_to_matrix_r2(weight_table)
//
// # Using `bench` package to evaluate the two versions of this function
//   if (!requireNamespace("bench", quietly = TRUE)) {
//     install.packages("bench")
//   }
//   bench::mark(table_to_matrix_r(weight_table),
//               table_to_matrix(weight_table)[genes, genes],
//               table_to_matrix_r2(weight_table))
