#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix DT2Matrix(DataFrame weightDT) {
  CharacterVector regulator = weightDT[0];
  CharacterVector target = weightDT[1];
  NumericVector weight = weightDT[2];

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

// #' R code
// #' @title DT2MatrixR
// #'  Switch weight data table to weight matrix
// #' @param weightDT The weight data table of network
// #'
// #' @return Weight matrix
// #' @export
// #'
// DT2MatrixR <- function(weightDT) {
//   colnames(weightDT) <- c("regulator", "target", "weight")
//   genes <- gtools::mixedsort(unique(c(weightDT$regulator, weightDT$target)))
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
//
// DT2MatrixR2 <- function(weightDT) {
//   weightMatrix <- reshape2::acast(weightDT, regulator ~ target)[genes, genes]
//   weightMatrix[is.na(weightMatrix)] <- 0
//   return(weightMatrix)
// }
//
// #' @examples
// Rcpp::sourceCpp("src/DT2Matrix.cpp")
//   library(inferCSN)
//   data("exampleMatrix")
//   weightDT <- inferCSN(exampleMatrix, verbose = TRUE)
//   weightMatrix <- DT2Matrix(weightDT)
//   genes <- gtools::mixedsort(unique(c(weightDT$regulator, weightDT$target)))
//   weightMatrix <- weightMatrix[genes, genes]
// weightMatrixR <- DT2MatrixR(weightDT)
//   weightMatrixR2 <- DT2MatrixR2(weightDT)
//
// # Using `bench` package to evaluate the two versions of this function
//   bench::mark(DT2MatrixR(weightDT),
//               DT2Matrix(weightDT)[genes, genes],
//                                  DT2MatrixR2(weightDT))
