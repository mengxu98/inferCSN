#' compute.gene.rank
#' @details Function to compute page rank of TF+target networks
#'
#' @param weightList Result of GRN reconstruction
#' @param directedGraph If GRN is directed or not
#'
#' @return A data.table with three columns.
#' @export
#'
#' @examples
#' data("exampleDataMatrix")
#' weightList <- inferCSN(exampleDataMatrix)
#' compute.gene.rank(weightList)
compute.gene.rank <- function(weightList,
                              directedGraph = FALSE) {
  if (!is.null(weightList)) {
    if (nrow(weightList) > 3) weightList <- weightList[, 1:3]
    colnames(weightList) <- c("regulatoryGene", "targetGene", "weight")
  } else {
    stop("Please input data......")
  }
  tfnet <- igraph::graph_from_data_frame(weightList, directed = directedGraph)
  pageRank <- data.frame(igraph::page_rank(tfnet, directed = directedGraph)$vector)
  colnames(pageRank) <- c("pageRank")
  pageRank$gene <- rownames(pageRank)
  pageRank <- pageRank[, c("gene", "pageRank")]
  pageRank <- pageRank[order(pageRank$pageRank, decreasing = TRUE), ]
  pageRank$is_regulator <- FALSE
  pageRank$is_regulator[pageRank$gene %in% unique(weightList$regulatoryGene)] <- TRUE
  return(pageRank)
}
