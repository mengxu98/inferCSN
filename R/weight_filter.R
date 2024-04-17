#' @title weight_filter
#'
#' @param weight_table weight_table
#' @param method method
#'
#' @return Filtered weight table
#' @export
#'
#' @examples
#' library(inferCSN)
#' data("example_matrix")
#' data("example_ground_truth")
#' weight_table <- inferCSN(example_matrix, verbose = TRUE)
#' weight_table_new <- weight_filter(weight_table)
#' network.heatmap(
#'   example_ground_truth[, 1:3],
#'   heatmap_title = "Ground truth",
#'   show_names = TRUE,
#'   rect_color = "gray90"
#' )
#' network.heatmap(
#'   weight_table,
#'   heatmap_title = "Raw",
#'   show_names = TRUE,
#'   rect_color = "gray90"
#' )
#' network.heatmap(
#'   weight_table_new,
#'   heatmap_title = "Filtered",
#'   show_names = TRUE,
#'   rect_color = "gray90"
#' )
#'
#' auc.calculate(
#'   weight_table,
#'   example_ground_truth,
#'   plot = TRUE
#'  )
#' auc.calculate(
#'   weight_table_new,
#'   example_ground_truth,
#'   plot = TRUE
#'  )
weight_filter <- function(
    weight_table,
    method = "max") {
  weight_table <- weight_table
  weight_table$edge <- paste(
    weight_table$regulator,
    weight_table$target,
    sep = "_"
  )

  weight_table_new <- data.frame(
    edge = paste(weight_table$target, weight_table$regulator, sep = "_"),
    weight = weight_table$weight
  )
  rownames(weight_table_new) <- weight_table_new$edge
  weight_table_new <- weight_table_new[weight_table$edge, ]
  weight_table$weight_new <- weight_table_new$weight
  weight_table <- dplyr::filter(weight_table, abs(weight) > abs(weight_new))
  weight_table <- weight_table[, c(1:3)]

  return(weight_table)
}
