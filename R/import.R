#' @import ggplot2 ggraph ggnetwork
#'
#' @importFrom Rcpp evalCpp sourceCpp
#' @importFrom RcppArmadillo armadillo_version
#' @importFrom RcppParallel RcppParallelLibs
#' @importFrom stats coef predict
NULL

utils::globalVariables(
  c(
    "Actual",
    "Category",
    "celltype",
    "cluster",
    "count",
    "curvetype",
    "degree",
    "edges",
    "from",
    "i",
    "id",
    "Interaction",
    "label_genes",
    "Metric",
    "name",
    "P_value",
    "P_value_contrary",
    "Predicted",
    "regulator",
    "targets_num",
    "target",
    "to",
    "type",
    "Value",
    "weight",
    "x",
    "xend",
    "y",
    "yend"
  )
)
