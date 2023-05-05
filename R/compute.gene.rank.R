#' @title compute.gene.rank
#' @details Function to compute page rank of TF+target networks
#'
#' @param weightDT The weight data table of network.
#' @param directedGraph If GRN is directed or not
#'
#' @return A data.table with three columns
#' @export
#'
#' @examples
#' data("exampleDataMatrix")
#' weightDT <- inferCSN(exampleDataMatrix)
#' ranks <- compute.gene.rank(weightDT)
#' head(ranks)
compute.gene.rank <- function(weightDT,
                              directedGraph = FALSE) {
  if (!is.null(weightDT)) {
    if (nrow(weightDT) > 3) weightDT <- weightDT[, 1:3]
    colnames(weightDT) <- c("regulatoryGene", "targetGene", "weight")
  } else {
    stop("Please input data......")
  }
  tfnet <- igraph::graph_from_data_frame(weightDT, directed = directedGraph)
  pageRank <- data.frame(igraph::page_rank(tfnet, directed = directedGraph)$vector)
  colnames(pageRank) <- c("pageRank")
  pageRank$gene <- rownames(pageRank)
  pageRank <- pageRank[, c("gene", "pageRank")]
  pageRank <- pageRank[order(pageRank$pageRank, decreasing = TRUE), ]
  pageRank$is_regulator <- FALSE
  pageRank$is_regulator[pageRank$gene %in% unique(weightDT$regulatoryGene)] <- TRUE
  return(pageRank)
}
