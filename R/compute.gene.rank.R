#' @title Compute and rank TFs in network
#'
#' @param weightDT The weight data table of network
#' @param directedGraph If GRN is directed or not
#'
#' @return A data.table with three columns
#' @export
#'
#' @examples
#' library(inferCSN)
#' data("exampleMatrix")
#' weightDT <- inferCSN(exampleMatrix)
#' ranks <- compute.gene.rank(weightDT)
#' head(ranks)
#'
compute.gene.rank <- function(weightDT,
                              directedGraph = FALSE) {
  colnames(weightDT) <- c("regulatory", "target", "weight")
  weightDT$weight <- abs(weightDT$weight)
  tfnet <- igraph::graph_from_data_frame(weightDT, directed = directedGraph)
  pageRank <- data.frame(igraph::page_rank(tfnet, directed = directedGraph)$vector)
  colnames(pageRank) <- c("pageRank")
  pageRank$gene <- rownames(pageRank)
  pageRank <- pageRank[, c("gene", "pageRank")]
  pageRank <- pageRank[order(pageRank$pageRank, decreasing = TRUE), ]
  pageRank$isRegulator <- FALSE
  pageRank$isRegulator[pageRank$gene %in% unique(weightDT$regulatory)] <- TRUE
  return(pageRank)
}
