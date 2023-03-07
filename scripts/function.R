#' loadpackages
#'  load packages
#' @param ... packages ...length()
#'
#' @return NULL
#' @export
#'
#' @examples load.packages(Seurat,ggplot2)
load.packages = function(...) {
  print(...length())
  pkgs = as.character(substitute(list(...)))[-1]   # 将输入转化为字符向量
  # print(pkgs)
  for (pkg in pkgs) {
    # 如果库里没有这个包,则去下载
    if (!requireNamespace(pkg, quietly = TRUE)) {
      # 加载 (传参%>%)、(网页爬取) 的包
      if (!requireNamespace("dplyr", quietly = TRUE)) {
        install.packages("dplyr")
      }
      library(dplyr)
      if (!requireNamespace("rvest", quietly = TRUE)) {
        install.packages("rvest")
      }
      library(rvest)
      # 爬取cran和bioconductor中包的列表
      message("[", Sys.time(), "] -----: No package: ", pkg, " in R environment!")
      CRANpackages = available.packages() %>%
        as.data.frame() %>%
        select(Package) %>%
        mutate(source = "CRAN")
      url = "https://www.bioconductor.org/packages/release/bioc/"
      biocPackages = url %>%
        read_html() %>%
        html_table() %>%
        .[[1]] %>%
        select(Package) %>%
        mutate(source = "BioConductor")
      # 如果在cran库中
      if (pkg %in% CRANpackages$Package) {
        message("[", Sys.time(), "] -----: Now install package: ", pkg, " from CRAN!")
        install.packages(pkg)
        library(pkg, character.only = T)
      }
      # 如果在bioconductor库中
      else if (pkg %in% biocPackages$Package) {
        message("[", Sys.time(), "] -----: Now install package: ", pkg, " from BioConductor!")
        BiocManager::install(pkg)
        library(pkg, character.only = T)
      }
      # 都不在,去github中查找
      else {
        # 加载包
        if (!requireNamespace("githubinstall", quietly = TRUE)) {
          install.packages("githubinstall")
        }
        library(githubinstall)
        # githubinstall(package)
        # 查找这个包的信息并输出
        gh_suggest(pkg)
      }
    }
    # 库里有这个包,直接加载
    else {
      message("[", Sys.time(), "] -----: ", pkg, " is in R environment!")
      library(pkg, character.only = T)
    }
  }
}

#' run.seurat.pipeline
#'
#' @param seurat_object
#' @param nfeatures
#' @param dims_Neighbors
#' @param dims_UMAP
#' @param dims_TSNE
#' @param resolution
#' @param verbose
#'
#' @return
#' @export
#'
#' @examples
run.seurat.pipeline = function(seurat_object,
                               nfeatures = 2000,
                               dims_Neighbors = 1:10,
                               dims_UMAP = 1:10,
                               dims_TSNE = 1:10,
                               resolution = 1.0,
                               verbose = FALSE){
  message("Running Seurat!")
  library(Seurat)
  seurat_object = NormalizeData(seurat_object, verbose = verbose)
  seurat_object = FindVariableFeatures(seurat_object, nfeatures = nfeatures, verbose = verbose)
  seurat_object = ScaleData(seurat_object, features = VariableFeatures(object = seurat_object), verbose = verbose)
  seurat_object = RunPCA(seurat_object, features = VariableFeatures(object = seurat_object), verbose = verbose)
  seurat_object = FindNeighbors(seurat_object, dims = dims_Neighbors, verbose = verbose)
  seurat_object = FindClusters(seurat_object, resolution = resolution, verbose = verbose)
  seurat_object = RunUMAP(seurat_object, dims = dims_UMAP, verbose = verbose)
  seurat_object = RunTSNE(seurat_object, dims = dims_TSNE, verbose = verbose)
  message("Finish Seurat!")
  return(seuratObject)
}


#### 将一个大的稀疏矩阵转换成矩阵  ####
Rcpp::sourceCpp(code='
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
IntegerMatrix asMatrix(NumericVector rp, NumericVector cp, NumericVector z, int nrows, int ncols)
{
  int k = z.size();
  IntegerMatrix  mat(nrows, ncols);
  for (int i = 0; i < k; i++)
  {
    mat(rp[i],cp[i]) = z[i];
  }
  return mat;
}
')
as_matrix = function(matrix) {
  row_pos = matrix@i
  col_pos = findInterval(seq(matrix@x)-1,matrix@p[-1])
  tmp = asMatrix(rp = row_pos, cp = col_pos, z = matrix@x,
                 nrows =  matrix@Dim[1], ncols = matrix@Dim[2])
  row.names(tmp) = matrix@Dimnames[[1]]
  colnames(tmp) = matrix@Dimnames[[2]]
  return(tmp)
}


#' annotation.celltype
#'
#' @param seu_obj
#' @param method
#'
#' @return
#' @export
#'
#' @examples
annotation.celltype = function(seuObj, method = "celltypist") {
  if (method == "celltypist") {
    loadpackages(reticulate)
    use_condaenv(condaenv = "scRNA", required = TRUE, conda = "/usr/local/anaconda3/bin/conda")
    pandas = import("pandas")
    numpy = import("numpy")
    scanpy = import("scanpy")
    celltypist = import("celltypist")
    tryCatch(
      {
        raw.matrix = t(as.matrix(seu_obj[["RNA"]]@counts))
      },
      warning = function(w) {
        print("warning:矩阵太大,使用as_matrix方法")
        raw.matrix = t(as_matrix(seu_obj[["RNA"]]@counts))
      },
      error = function(e) {
        print("error:矩阵太大,使用as_matrix方法")
        raw.matrix = t(as_matrix(seu_obj[["RNA"]]@counts))
      }
    )
    message("[", Sys.time(), "] -----: Run 'celltypist'!")
    adata = scanpy$AnnData(
      X = numpy$array(raw.matrix),
      obs = pandas$DataFrame(seu_obj@meta.data),
      var = pandas$DataFrame(data.frame(
        gene = rownames(seu_obj[["RNA"]]@counts),
        row.names = rownames(seu_obj[["RNA"]]@counts)
      ))
    )
    scanpy$pp$normalize_total(adata, target_sum = 1e4)
    scanpy$pp$log1p(adata)
    predictions = celltypist$annotate(adata, model = "Cells_Lung_Airway.pkl", majority_voting = FALSE)
    seu_obj$celltypist.labels = predictions$predicted_labels
  }
  else if (method == "SingleR") {
    # 加载包
    loadpackages(SingleR,Seurat)
    # 下载参考数据集
    # ref.data <- HumanPrimaryCellAtlasData(ensembl=TRUE)
    # 加载本地的参考数据集
    load("/media/ssd/cx-workplace/data/SingleR_labels/HumanPrimaryCellAtlas_hpca.se_human.RData")
    normalized.matrix = GetAssayData(seu_obj, slot="data")  # 获取归一化表达矩阵
    #clusters.labels = seu_obj$seurat_clusters              # 获取聚类信息(似乎已弃用)
    message("[", Sys.time(), "] -----: Run 'SingleR'!")
    predicted_label = SingleR(test = normalized.matrix,
                              ref = hpca.se,
                              # method = "cluster",
                              # clusters = clusters.labels,
                              labels = hpca.se$label.main) # 执行singleR主程序
    seu_obj$SingleR.labels = predicted_label$labels        # 添加到meta.data
  }
  return(seu_obj)
}

