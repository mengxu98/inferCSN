#' @title data.processing
#'
#' @param data [Default = NULL] A matrix, data table, Seurat or SingleCellExperiment object
#' @param normalize [Default = FALSE] Data normalize
#' @param verbose [Default = FALSE] Print detailed information
#'
#' @return Matrix
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

  # Remove ribo.genes
  genes.all <- rownames(matrix)
  ribo.genes <- genes.all[grep('^RPS[0-9]*|^RPL[0-9]*', genes.all)]
  mito.genes <- genes.all[grep('^MRPS[0-9]*|^MRPL[0-9]*', genes.all)]
  nduf.genes <- genes.all[grep('^NDUF', genes.all)]
  matrix <- matrix[setdiff(genes.all, c(ribo.genes, mito.genes, nduf.genes)),]

  if (normalize) {

    # Norm data
    geneFiltered <- apply(matrix, 1, function(x) {
      sum(x >= 1) >= 1
    })

    Sizes <- EBSeq::MedianNorm(matrix)
    if (is.na(Sizes[1])) {
      matrixMedNorm <- matrix
      # message("alternative normalization method is applied")
      # Sizes <- EBSeq::MedianNorm(matrix, alternative=TRUE)
    } else {
      matrixMedNorm <- EBSeq::GetNormalizedMat(matrix, Sizes)
    }

    matrix <- log2(matrixMedNorm + 1)
    matrix <- matrix[which(geneFiltered), ]
  }
  matrix <- as.data.frame(t(matrix))
  return(matrix)
}

#' @title smooths expression
#' @description
#'  smooths expression across metadata in path
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
