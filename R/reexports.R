#' @import methods
#'
#' @importClassesFrom Matrix sparseMatrix
#' @importFrom Rcpp evalCpp sourceCpp
#' @importFrom stats coef predict
#' @importFrom utils head tail
NULL

utils::globalVariables(
  c(
    "regulator",
    "target"
  )
)
