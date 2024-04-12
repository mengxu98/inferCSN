#' @title Calculate and rank TFs in network
#'
#' @inheritParams net.format
#' @param directed If network is directed or not.
#'
#' @return A data.table with three columns
#' @export
#'
#' @examples
#' library(inferCSN)
#' data("example_matrix")
#' weight_table <- inferCSN(example_matrix)
#' head(calculate.gene.rank(weight_table))
#' head(calculate.gene.rank(weight_table, regulators = "g1"))
calculate.gene.rank <- function(
    weight_table,
    regulators = NULL,
    targets = NULL,
    directed = FALSE) {
  colnames(weight_table) <- c("regulator", "target", "weight")
  weight_table <- net.format(
    weight_table,
    regulators = regulators,
    targets = targets
  )

  network <- igraph::graph_from_data_frame(
    weight_table,
    directed = directed
  )
  page_rank_res <- data.frame(
    igraph::page_rank(network, directed = directed)$vector
  )
  colnames(page_rank_res) <- c("page_rank")
  page_rank_res$gene <- rownames(page_rank_res)
  page_rank_res <- page_rank_res[, c("gene", "page_rank")]
  page_rank_res <- page_rank_res[order(
    page_rank_res$page_rank,
    decreasing = TRUE
  ), ]
  page_rank_res$is_regulator <- FALSE
  page_rank_res$is_regulator[
    page_rank_res$gene %in% unique(
      weight_table$regulator
    )
  ] <- TRUE

  return(page_rank_res)
}
