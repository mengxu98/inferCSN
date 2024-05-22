#' @import Matrix ggplot2 patchwork
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
    "i",
    "Interaction",
    "name",
    "degree",
    "edges",
    "curvetype",
    "weight_new",
    "celltype",
    "from",
    "id",
    "label_genes",
    "targets_num",
    "to"
  )
)
