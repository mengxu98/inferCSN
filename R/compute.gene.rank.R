#' compute.gene.rank
#' @details Function to compute page rank of TF+target networks
#'
#' @param weightdf Result of GRN reconstruction
#' @param directedGraph If GRN is directed or not
#'
#' @return
#' @export
#'
#' @examples
compute.gene.rank <- function(weightdf, directedGraph = FALSE) {
  if (!is.null(weightdf)) {
    if (nrow(weightdf)==3) {
      colnames(weightdf) <- c("regulatoryGene", "targetGene", "weight")
    }else{
      weightdf <- weightdf[, 1:3]
      colnames(weightdf) <- c("regulatoryGene", "targetGene", "weight")
    }
  }
  tfnet <- igraph::graph_from_data_frame(weightdf, directed = directedGraph)
  pageRank <- data.frame(igraph::page_rank(tfnet, directed = directedGraph)$vector)
  colnames(pageRank) <- c("pageRank")
  pageRank$gene <- rownames(pageRank)
  pageRank <- pageRank[, c("gene", "pageRank")]
  pageRank <- pageRank[order(pageRank$pageRank, decreasing = TRUE), ]
  pageRank$is_regulator <- FALSE
  pageRank$is_regulator[pageRank$gene %in% unique(weightdf$regulatoryGene)] <- TRUE
  return(pageRank)
}
