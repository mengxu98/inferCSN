#' @title dynamic.genes
#' @description
#'  Find genes expressed dynamically using slingshot approach
#'
#' @param matrix properly normalized expression matrix
#' @param metadata sample table that includes pseudotime, rownames = cells, and a group column
#' @param path vector of group names to include
#' @param method method to find dynamic genes. Defaults to "gam", if "tradeseq", must provide matrixRaw
#' @param matrixRaw raw expression matrix, required if method="tradeseq"
#'
#' @import tradeSeq
#' @return pvals and cell info
#'
#' @export
#'
dynamic.genes <- function(matrix,
                          metadata,
                          path = NULL,
                          method = "gam",
                          matrixRaw = NULL) {
  if (method == "gam") {
    pseudotime <- metadata$pseudotime
    names(pseudotime) <- as.vector(metadata$cells)
    cat("starting gammma...\n")
    gpChr <- gam.fit(matrix[, names(pseudotime)], rownames(matrix), pseudotime)
    gpChr <- gam.fit(matrix, rownames(matrix), pseudotime)
    cells <- data.frame(cells = names(pseudotime), pseudotime = pseudotime)
    rownames(cells) <- names(pseudotime)
    cells <- cells[order(cells$pseudotime), ]
    ans <- list(genes = gpChr, cells = cells)
  } else if (method == "tradeseq"){
    if (is.null(matrixRaw)) {
      matrixRaw <- matrix
      # stop("Must provide matrixRaw for TradeSeq.")
    }
    # subset raw data based on normalized data
    matrixRaw <- matrixRaw[, rownames(metadata)]
    matrixRaw <- matrixRaw[rownames(matrix), ]
    pt <- as.data.frame(metadata[, "pseudotime"])
    rownames(pt) <- rownames(metadata)
    colnames(pt) <- "pseudotime"
    cw <- as.matrix(rep(1, nrow(pt)))
    rownames(cw) <- rownames(metadata)
    ts <- tradeSeq::fitGAM(
      as.matrix(matrixRaw),
      pseudotime = as.matrix(pt),
      cellWeights = cw
    )
    ATres <- associationTest(ts) # %>% na.omit()
    genes <- ATres$pvalue
    names(genes) <- rownames(ATres)
    cells <- data.frame(
      cells = metadata$cells,
      pseudotime = metadata$pseudotime
    )
    rownames(cells) <- metadata$cells
    cells <- cells[order(cells$pseudotime), ]
    ans <- list(genes = genes, cells = cells)
  }
  ans
}

#' @title gam.fit
#'
#' @param matrix matrix
#' @param genes genes
#' @param celltime celltime
#'
#' @return list
#' @export
#'
gam.fit <- function(matrix,
                    genes,
                    celltime) {
  genesInter <- intersect(genes, rownames(matrix))
  # could print out if any missing
  ans <- apply(matrix[genesInter, names(celltime)], 1, function(z) {
    d <- data.frame(z = z, t = celltime)
    tmp <- gam::gam(z ~ gam::lo(celltime), data = d)
    p <- summary(tmp)[4][[1]][1, 5]
    p
  })
  ans
}
