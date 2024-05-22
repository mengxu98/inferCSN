#' @title Calculate and rank TFs in network
#'
#' @inheritParams network_format
#' @param directed If network is directed or not.
#'
#' @return A data.table with three columns
#' @export
#'
#' @examples
#' data("example_matrix")
#' network_table <- inferCSN(example_matrix)
#' head(calculate.gene.rank(network_table))
#' head(calculate.gene.rank(network_table, regulators = "g1"))
calculate.gene.rank <- function(
    network_table,
    regulators = NULL,
    targets = NULL,
    directed = FALSE) {
  colnames(network_table) <- c("regulator", "target", "weight")
  network_table <- network_format(
    network_table,
    regulators = regulators,
    targets = targets
  )

  network <- igraph::graph_from_data_frame(
    network_table,
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
      network_table$regulator
    )
  ] <- TRUE

  return(page_rank_res)
}
