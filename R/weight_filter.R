#' @title weight_filter
#'
#' @param network_table network_table
#' @param method method
#'
#' @return Filtered weight table
#' @export
#'
#' @examples
#' data("example_matrix")
#' data("example_ground_truth")
#' network_table <- inferCSN(example_matrix)
#' network_table_new <- weight_filter(network_table)
#' network.heatmap(
#'   example_ground_truth[, 1:3],
#'   heatmap_title = "Ground truth",
#'   show_names = TRUE,
#'   rect_color = "gray90"
#' )
#' network.heatmap(
#'   network_table,
#'   heatmap_title = "Raw",
#'   show_names = TRUE,
#'   rect_color = "gray90"
#' )
#' network.heatmap(
#'   network_table_new,
#'   heatmap_title = "Filtered",
#'   show_names = TRUE,
#'   rect_color = "gray90"
#' )
#'
#' auc.calculate(
#'   network_table,
#'   example_ground_truth,
#'   plot = TRUE
#' )
#' auc.calculate(
#'   network_table_new,
#'   example_ground_truth,
#'   plot = TRUE
#' )
weight_filter <- function(
    network_table,
    method = "max") {
  network_table$edge <- paste(
    network_table$regulator,
    network_table$target,
    sep = "_"
  )

  network_table_new <- data.frame(
    edge = paste(network_table$target, network_table$regulator, sep = "_"),
    weight = network_table$weight
  )
  rownames(network_table_new) <- network_table_new$edge
  network_table_new <- network_table_new[network_table$edge, ]
  network_table$weight_new <- network_table_new$weight
  if (method == "max") {
    network_table <- dplyr::filter(network_table, abs(weight) > abs(weight_new))
  }
  network_table <- network_table[, c("regulator", "target", "weight")]

  return(network_table)
}
