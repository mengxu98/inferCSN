#' data.processing
#'
#' @param data An object
#' @param normalize Data normalize
#' @param verbose Print detailed information
#'
#' @return Matrix
#' @export
#'
data.processing <- function(data,
                            normalize = normalize,
                            verbose = verbose){
  if (class(data)[1] == "matrix") {
    matrix <- t(data)
  } else if (class(data)[1] == "data.frame") {
    matrix <- t(data)
  } else if (class(data)[1] == "Seurat") {
    matrix <- t(as.matrix(data@assays$RNA@counts))
  } else if (class(data)[1] == "SingleCellExperiment") {
    matrix <- t(as.matrix(sce@assays@data$counts))
  } else {
    stop("Error")
  }

  if (verbose) message("Data processing......")
  # Remove ribo.genes
  genes.all <- rownames(matrix)
  ribo.genes <- genes.all[grep('^RPS[0-9]*|^RPL[0-9]*', genes.all)]
  mito.genes <- genes.all[grep('^MRPS[0-9]*|^MRPL[0-9]*', genes.all)]
  nduf.genes <- genes.all[grep('^NDUF', genes.all)]
  matrix <- matrix[setdiff(genes.all, c(ribo.genes, mito.genes, nduf.genes)), ]

  if (normalize) {

    # Norm data
    geneFiltered <- apply(matrix,1,function(x){
      sum(x >= 1) >= 1
    })

    Sizes <- EBSeq::MedianNorm(matrix)
    if(is.na(Sizes[1])){
      matrixMedNorm <- matrix
      # message("alternative normalization method is applied")
      # Sizes <- EBSeq::MedianNorm(matrix, alternative=TRUE)
    }else{
      matrixMedNorm <- EBSeq::GetNormalizedMat(matrix, Sizes)
    }

    matrix <- log2(matrixMedNorm + 1)
    matrix <- matrix[which(geneFiltered),]
  }
  matrix <- as.data.frame(t(matrix))
  return(matrix)
}
