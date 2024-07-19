#' @title Rank TFs and genes in network
#'
#' @inheritParams network_format
#' @param directed Whether the network is directed.
#'
#' @return A table of gene rank.
#' @export
#'
#' @examples
#' data("example_matrix")
#' network_table <- inferCSN(example_matrix)
#' head(calculate_gene_rank(network_table))
#' head(calculate_gene_rank(network_table, regulators = "g1"))
#' head(calculate_gene_rank(network_table, targets = "g1"))
calculate_gene_rank <- function(
    network_table,
    regulators = NULL,
    targets = NULL,
    directed = FALSE) {
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
  colnames(page_rank_res) <- c("rank_value")
  page_rank_res$gene <- rownames(page_rank_res)
  page_rank_res <- page_rank_res[, c("gene", "rank_value")]
  page_rank_res <- page_rank_res[order(
    page_rank_res$rank_value,
    decreasing = TRUE
  ), ]
  page_rank_res$regulator <- ifelse(
    page_rank_res$gene %in% unique(network_table$regulator),
    "TRUE", "FALSE"
  )

  return(page_rank_res)
}
