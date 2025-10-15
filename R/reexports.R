#' @import ggplot2 ggraph ggnetwork
#'
#' @importClassesFrom Matrix sparseMatrix
#' @importFrom Rcpp evalCpp sourceCpp
#' @importFrom stats coef predict
#' @importFrom thisutils `%ss%`
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
