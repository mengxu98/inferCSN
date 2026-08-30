#' @title Subsample an expression matrix
#'
#' @param matrix Input matrix.
#' @param subsampling_method Subsampling strategy.
#' @param subsampling_ratio Fraction retained or aggregated.
#' @param seed Random seed.
#' @param verbose Whether to report progress.
#' @param ... Additional arguments.
#' @return The subsampled matrix.
#' @export
subsampling <- function(
  matrix,
  subsampling_method = c("sample", "meta_cells", "pseudobulk"),
  subsampling_ratio = 1,
  seed = 1,
  verbose = TRUE,
  ...
) {
  subsampling_method <- match.arg(
    subsampling_method
  )

  if (!(is.numeric(subsampling_ratio) && subsampling_ratio > 0 && subsampling_ratio <= 1)) {
    thisutils::log_message(
      "Please set {.arg subsampling_ratio} value between: (0, 1]",
      message_type = "error"
    )
  }
  if (subsampling_ratio >= 1) {
    return(matrix)
  }
  if (methods::is(matrix, "sparseMatrix")) {
    return_sparse <- TRUE
  } else {
    return_sparse <- FALSE
  }

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
      pseudobulk(
        matrix = matrix,
        ratio = subsampling_ratio,
        ...
      )
    }
  )

  if (return_sparse) {
    matrix <- thisutils::as_matrix(matrix, return_sparse = TRUE)
  } else {
    matrix <- thisutils::as_matrix(matrix)
  }

  thisutils::log_message(
    "Subsample matrix generated, dimensions: {.val {nrow(matrix)}} cells by {.val {ncol(matrix)}} genes",
    verbose = verbose
  )

  return(matrix)
}

pseudobulk <- function(
  matrix,
  ratio = 0.5,
  k = 50,
  seed = 1,
  prefix = "pseudobulk_",
  ...
) {
  n_samples <- round(nrow(matrix) * ratio)

  knn_res <- build_knn(
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
  rownames(agg_matrix) <- paste0(prefix, seq_len(n_samples))

  return(agg_matrix)
}
