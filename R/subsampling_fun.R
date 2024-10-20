#' @title Subsampling function
#'
#' @description
#' This function subsamples a matrix using either random sampling or meta cells method.
#'
#' @inheritParams inferCSN
#' @param matrix The input matrix to be subsampled.
#'
#' @return The subsampled matrix.
#'
#' @export
#'
#' @examples
#' data("example_matrix")
#' data("example_ground_truth")
#' subsample_matrix <- subsampling_fun(
#'   example_matrix,
#'   subsampling_ratio = 0.5
#' )
#' subsample_matrix_2 <- subsampling_fun(
#'   example_matrix,
#'   subsampling_method = "meta_cells",
#'   subsampling_ratio = 0.5,
#'   fast_pca = FALSE
#' )
#'
#' calculate_auc(
#'   inferCSN(example_matrix),
#'   example_ground_truth,
#'   plot = TRUE
#' )
#' calculate_auc(
#'   inferCSN(subsample_matrix),
#'   example_ground_truth,
#'   plot = TRUE
#' )
#' calculate_auc(
#'   inferCSN(subsample_matrix_2),
#'   example_ground_truth,
#'   plot = TRUE
#' )
subsampling_fun <- function(
    matrix,
    subsampling_method = "sample",
    subsampling_ratio = 1,
    seed = 1,
    verbose = TRUE,
    ...) {
  if (subsampling_ratio >= 1) {
    return(matrix)
  }
  if (methods::is(matrix, "sparseMatrix")) {
    return_sparse <- TRUE
  } else {
    return_sparse <- FALSE
  }

  subsampling_method <- match.arg(
    subsampling_method,
    c("sample", "meta_cells")
  )

  set.seed(seed)
  matrix <- switch(
    EXPR = subsampling_method,
    "sample" = {
      sample_count <- nrow(matrix)
      subsample_count <- round(sample_count * subsampling_ratio)
      matrix[sample(sample_count, subsample_count), ]
    },
    "meta_cells" = {
      meta_cells(
        matrix = matrix,
        gamma = 1 / subsampling_ratio,
        ...
      )
    }
  )

  if (return_sparse) {
    matrix <- as_matrix(matrix, sparse = TRUE)
  } else {
    matrix <- as_matrix(matrix)
  }

  log_message(
    "Subsample matrix generated, dimensions: ",
    nrow(matrix), " cells by ",
    ncol(matrix), " genes.",
    verbose = verbose
  )

  return(matrix)
}
