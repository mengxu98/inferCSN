#' @import Matrix ggplot2 ggraph patchwork
#'
#' @importFrom Rcpp evalCpp sourceCpp
#' @importFrom RcppArmadillo armadillo_version
#' @importFrom RcppParallel RcppParallelLibs
#' @importFrom stats coef predict
#' @importFrom utils methods
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
    "weight_new",
    "celltype",
    "from",
    "id",
    "label_genes",
    "targets_num",
    "to",
    "P_value",
    "P_value_contrary",
    "v",
    "v_new"
  )
)
