#' @import Matrix ggplot2 patchwork
#'
#' @importFrom Rcpp evalCpp
#' @importFrom RcppArmadillo armadillo_version
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
