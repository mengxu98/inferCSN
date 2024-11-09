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
#' subsample_matrix <- subsampling(
#'   example_matrix,
#'   subsampling_ratio = 0.5
#' )
#' subsample_matrix_2 <- subsampling(
#'   example_matrix,
#'   subsampling_method = "meta_cells",
#'   subsampling_ratio = 0.5,
#'   fast_pca = FALSE
#' )
#' subsample_matrix_3 <- subsampling(
#'   example_matrix,
#'   subsampling_method = "pseudobulk",
#'   subsampling_ratio = 0.5
#' )
#'
#' calculate_metrics(
#'   inferCSN(example_matrix),
#'   example_ground_truth,
#'   plot = TRUE
#' )
#' calculate_metrics(
#'   inferCSN(subsample_matrix),
#'   example_ground_truth,
#'   plot = TRUE
#' )
#' calculate_metrics(
#'   inferCSN(subsample_matrix_2),
#'   example_ground_truth,
#'   plot = TRUE
#' )
#' calculate_metrics(
#'   inferCSN(subsample_matrix_3),
#'   example_ground_truth,
#'   plot = TRUE
#' )
subsampling <- function(
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
    c("sample", "meta_cells", "pseudobulk")
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
    },
    "pseudobulk" = {
      .pseudobulk(
        matrix = matrix,
        ratio = subsampling_ratio,
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

.pseudobulk <- function(
    matrix,
    ratio = 0.5,
    k = 50,
    seed = 1,
    ...) {
  n_samples <- round(nrow(matrix) * ratio)

  knn_res <- .build_knn(
    matrix = matrix,
    k = k,
    from = "coordinates",
    use_nn2 = TRUE,
    ...
  )

  set.seed(seed)
  seed_cells <- sample(seq_len(nrow(matrix)), n_samples)

  neighbors <- knn_res$idx

  agg_matrix <- matrix(
    0,
    nrow = n_samples,
    ncol = ncol(matrix)
  )
  for (i in seq_len(n_samples)) {
    cell_idx <- seed_cells[i]
    neighbor_idx <- neighbors[cell_idx, ]
    cells_to_aggregate <- c(cell_idx, neighbor_idx)
    cells_to_aggregate <- cells_to_aggregate[!is.na(cells_to_aggregate)]
    agg_matrix[i, ] <- colMeans(
      matrix[cells_to_aggregate, , drop = FALSE]
    )
  }

  colnames(agg_matrix) <- colnames(matrix)
  rownames(agg_matrix) <- paste0("pseudobulk_", seq_len(n_samples))

  return(agg_matrix)
}
