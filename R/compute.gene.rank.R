#' @title Compute and rank TFs in network
#'
#' @param weight_table The weight data table of network
#' @param directed If GRN is directed or not
#'
#' @return A data.table with three columns
#' @export
#'
#' @examples
#' library(inferCSN)
#' data("example_matrix")
#' weight_table <- inferCSN(example_matrix)
#' ranks <- compute.gene.rank(weight_table)
#' head(ranks)
compute.gene.rank <- function(
    weight_table,
    directed = FALSE) {
  colnames(weight_table) <- c("regulatory", "target", "weight")
  weight_table$weight <- abs(weight_table$weight)
  tfnet <- igraph::graph_from_data_frame(
    weight_table,
    directed = directed
  )
  page_rank_res <- data.frame(
    igraph::page_rank(tfnet, directed = directed)$vector
  )
  colnames(page_rank_res) <- c("pageRank")
  page_rank_res$gene <- rownames(page_rank_res)
  page_rank_res <- page_rank_res[, c("gene", "pageRank")]
  page_rank_res <- page_rank_res[order(
    page_rank_res$pageRank,
    decreasing = TRUE
  ), ]
  page_rank_res$isRegulator <- FALSE
  page_rank_res$isRegulator[
    page_rank_res$gene %in% unique(
      weight_table$regulatory
    )
  ] <- TRUE
  return(page_rank_res)
}
