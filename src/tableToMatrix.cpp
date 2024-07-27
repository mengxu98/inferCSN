#include <Rcpp.h>
#include <algorithm>

using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix tableToMatrix(DataFrame weight_table) {
  CharacterVector regulators = weight_table["regulator"];
  CharacterVector targets = weight_table["target"];
  NumericVector weight = weight_table["weight"];

  CharacterVector unique_regulators = unique(regulators);
  CharacterVector unique_targets = unique(targets);

  // // TODO:
  // // The two lines to sort genes, such as "g3", "g1", "g2" to "g1", "g2",
  // "g3",
  // // but these codes can not work successful, maybe return "g1.1", "g1.2",
  // "g3".
  // // Now, use R code: 'weight_matrix <- weight_matrix[unique_regulators,
  // unique_targets]',
  // // to temporarily overcome this problem
  // std::sort(unique_regulators.begin(), unique_regulators.end());
  // std::sort(unique_targets.begin(), unique_targets.end());

  int num_regulators = unique_regulators.size();
  int num_targets = unique_targets.size();

  NumericMatrix weight_matrix(num_regulators, num_targets);
  rownames(weight_matrix) = unique_regulators;
  colnames(weight_matrix) = unique_targets;

  for (int i = 0; i < weight_table.nrows(); ++i) {
    int regulators_index =
        std::distance(unique_regulators.begin(),
                      std::find(unique_regulators.begin(),
                                unique_regulators.end(), regulators[i]));
    int targets_index = std::distance(
        unique_targets.begin(),
        std::find(unique_targets.begin(), unique_targets.end(), targets[i]));
    weight_matrix(regulators_index, targets_index) = weight[i];
  }

  return weight_matrix;
}

/*

 # Older version, bug: can not handle data frame with unequal number regulators
 and targets. NumericMatrix tableToMatrix(DataFrame weight_table) {
 CharacterVector regulator = weight_table[0];
 CharacterVector target = weight_table[1];
 NumericVector weight = weight_table[2];

 CharacterVector genes = union_(regulator, target);

 int num_genes = genes.size();

 NumericMatrix weight_matrix(num_genes, num_genes);
 colnames(weight_matrix) = genes;
 rownames(weight_matrix) = genes;

 for (int i = 0; i < weight_table.nrows(); ++i) {
 int regulators_index = std::distance(genes.begin(), std::find(genes.begin(),
 genes.end(), regulator[i])); int targets_index = std::distance(genes.begin(),
 std::find(genes.begin(), genes.end(), target[i]));
 weight_matrix(regulators_index, targets_index) = weight[i];
 }

 return weight_matrix;
 }

 */

/*

#' R code version 1
#' @title table_to_matrix_v1
#'  Switch weight data table to weight matrix
#' @param weight_table The weight data table of network
#'
#' @return Weight matrix
#' @export
 table_to_matrix_v1 <- function(weight_table) {
 colnames(weight_table) <- c("regulator", "target", "weight")

 unique_regulators <- gtools::mixedsort(unique(weight_table$regulator))
 unique_targets <- gtools::mixedsort(unique(weight_table$target))

 weight_matrix <- matrix(
 0,
 nrow = length(unique_regulators),
 ncol = length(unique_targets),
 dimnames = list(unique_regulators, unique_targets)
 )

 for (i in 1:nrow(weight_table)) {
 weight_matrix[weight_table$regulator[i], weight_table$target[i]] <-
weight_table$weight[i]
 }

 return(weight_matrix)
 }

#' R code version 2
#' @title table_to_matrix_v2
#'  Switch weight data table to weight matrix
#' @param weight_table The weight data table of network
#'
#' @return Weight matrix
#' @export
 table_to_matrix_v2 <- function(weight_table) {
 unique_regulators <- gtools::mixedsort(unique(weight_table$regulator))
 unique_targets <- gtools::mixedsort(unique(weight_table$target))
 weight_matrix <- suppressMessages(
 reshape2::acast(
 weight_table,
 regulator ~ target,
 fill = 0
 )[unique_regulators, unique_targets]
 )

 return(weight_matrix)
 }

#' @examples
 Rcpp::sourceCpp("src/tableToMatrix.cpp")
 library(inferCSN)
 data("example_matrix")
 weight_table <- inferCSN(
 example_matrix,
 regulators = c("g1", "g2", "g3"),
 verbose = TRUE
 )

 unique_regulators <- gtools::mixedsort(unique(weight_table$regulator))
 unique_targets <- gtools::mixedsort(unique(weight_table$target))
 weight_matrix <- tableToMatrix(weight_table)
 weight_matrix <- weight_matrix[unique_regulators, unique_targets]
 weight_matrix_r1 <- table_to_matrix_v1(weight_table)
 weight_matrix_r2 <- table_to_matrix_v2(weight_table)

# Using `bench` package to evaluate different versions of this function
 if (!requireNamespace("bench", quietly = TRUE)) {
 install.packages("bench")
 }
 bench::mark(
 tableToMatrix(weight_table)[unique_regulators, unique_targets],
 table_to_matrix_v1(weight_table),
 table_to_matrix_v2(weight_table)
 )

 */
