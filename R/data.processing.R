#' @title data.processing
#'
#' @inheritParams inferCSN
#'
#' @return Normalized matrix
#' @export
#'
data.processing <- function(data,
                            normalize = FALSE,
                            verbose = FALSE) {
  if (class(data)[1] == "matrix") {
    matrix <- t(data)
  } else if (class(data)[1] == "data.frame") {
    matrix <- t(data)
  } else if (class(data)[1] == "Seurat") {
    matrix <- as.matrix(data@assays$RNA@counts)
  } else if (class(data)[1] == "SingleCellExperiment") {
    matrix <- as.matrix(sce@assays@data$counts)
  } else {
    stop("Error")
  }

  # Remove ribo.mito.nduf.genes
  genes.all <- rownames(matrix)
  ribo.genes <- genes.all[grep('^RPS[0-9]*|^RPL[0-9]*', genes.all)]
  mito.genes <- genes.all[grep('^MRPS[0-9]*|^MRPL[0-9]*', genes.all)]
  nduf.genes <- genes.all[grep('^NDUF', genes.all)]
  matrix <- matrix[setdiff(genes.all, c(ribo.genes, mito.genes, nduf.genes)),]

  # Normalize matrix
  if (normalize) {
    geneFiltered <- apply(matrix, 1, function(x) {
      sum(x >= 1) >= 1
    })
    sizes <- EBSeq::MedianNorm(matrix)
    if (is.na(sizes[1])) {
      matrixMedNorm <- matrix
    } else {
      matrixMedNorm <- EBSeq::GetNormalizedMat(matrix, sizes)
    }
    matrix <- log2(matrixMedNorm + 1)
    matrix <- matrix[which(geneFiltered), ]
  }
  matrix <- as.data.frame(t(matrix))
  return(matrix)
}

#' @title data.ksmooth
#' @description Smooths expression across pseudotime
#'
#' @param matrix  expression matrix
#' @param bandwidth bandwidth
#' @param metadata sample table that includes "cell"  "pseudotime" "group"
#'
#' @return smoothed matrix
#'
#' @export
#'
data.ksmooth <- function(matrix,
                         metadata,
                         bandwidth = 0.25) {
  bandwidth <- min(bandwidth, (max(metadata$pseudotime) - min(metadata$pseudotime)) / 10)

  pseudotime <- metadata$pseudotime
  names(pseudotime) <- as.vector(metadata$cells)
  sort(pseudotime)
  matrix <- matrix[, names(pseudotime)]
  matrixKS <- matrix(0, nrow = nrow(matrix), ncol = ncol(matrix))
  for (i in seq(nrow(matrix))) {
    yy <- ksmooth(pseudotime, matrix[i, ], kernel = "normal", bandwidth = bandwidth, x.points = pseudotime)
    matrixKS[i, ] <- yy$y
  }
  rownames(matrixKS) <- rownames(matrix)
  colnames(matrixKS) <- colnames(matrix)
  return(matrixKS)
}
