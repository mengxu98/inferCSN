#' @import ggplot2 ggraph ggnetwork
#'
#' @importFrom Rcpp evalCpp sourceCpp
#' @importFrom RcppArmadillo armadillo_version
#' @importFrom RcppParallel RcppParallelLibs
#' @importFrom stats coef predict
NULL

utils::globalVariables(
  c(
    "count",
    "x",
    "y",
    "xend",
    "yend",
    "regulator",
    "target",
    "weight",
    "i",
    "Interaction",
    "name",
    "degree",
    "edges",
    "cluster",
    "curvetype",
    "celltype",
    "from",
    "id",
    "label_genes",
    "targets_num",
    "to",
    "P_value",
    "P_value_contrary"
  )
)
