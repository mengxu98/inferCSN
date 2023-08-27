#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix DT2Matrix(DataFrame weightDT) {
  CharacterVector regulator = weightDT["regulator"];
  CharacterVector target = weightDT["target"];
  NumericVector weight = weightDT["weight"];

  CharacterVector genes = union_(regulator, target);

  int numGenes = genes.size();
  NumericMatrix weightMatrix(numGenes, numGenes);
  colnames(weightMatrix) = genes;
  rownames(weightMatrix) = genes;

  for (int i = 0; i < weightDT.nrows(); ++i) {
    int regIdx = std::distance(genes.begin(), std::find(genes.begin(), genes.end(), regulator[i]));
    int tarIdx = std::distance(genes.begin(), std::find(genes.begin(), genes.end(), target[i]));
    weightMatrix(regIdx, tarIdx) = weight[i];
  }

  return weightMatrix;
}

// R code
// #' @title Switch weight data table to weight matrix
// #'
// #' @param weightDT The weight data table of network
// #'
// #' @return Weight matrix
// #' @export
// #'
// DT2Matrix <- function(weightDT) {
//   colnames(weightDT) <- c("regulator", "target", "weight")
//   genes <- unique(c(weightDT$regulator, weightDT$target))
//
//   weightMatrix <- matrix(0,
//                          nrow = length(genes),
//                          ncol = length(genes),
//                          dimnames = list(genes, genes))
//
//   for (i in 1:nrow(weightDT)) {
//     weightMatrix[weightDT$regulator[i], weightDT$target[i]] <- weightDT$weight[i]
//   }
//
//   return(weightMatrix)
// }
