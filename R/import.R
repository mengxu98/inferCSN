#' @import Matrix
#'
#' @importFrom Rcpp evalCpp
#' @importFrom stats coef predict
#' @importFrom utils methods
NULL

utils::globalVariables(
  c(
    "x",
    "y",
    "xend",
    "yend",
    "regulator",
    "target",
    "weight",
    "Interaction",
    "name",
    "degree",
    "edges",
    "curvetype",
    "weight_new"
  )
)
