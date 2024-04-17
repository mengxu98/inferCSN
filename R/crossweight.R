#' Perform crossweighting
#'
#' @param weight_table GRN dataframe, the result of running reconstructargetRN or reconstructargetRN_GENIE3
#' @param matrix genes-by-cells expression matrix
#' @param meta_data result of running findDynGenes
#' @param lag lag window on which to run cross-correlation. Cross-correlaiton computed from -lag to +lag.
#' @param min minimum of weighting window. Edges with offsets (or absolute offsets if symmetric_filter=TRUE) less than min will not be negatively weighted.
#' @param max maximum of weighting window. Edges with offsets (or absolute offsets if symmetric_filter=TRUE) greater than max will have weights set to 0.
#' @param symmetric_filter whether or not to employ a symmetric weight scheme. If true, absolute offset is used in place of offset.
#' @param filter_thresh after crossweighting, edges with weights less than filter_thresh will be set to 0.
#'
#' @return weight_table with offset and weighted_score added
#'
#' @export
#' @examples
#' \dontrun{
#' library(inferCSN)
#' data("example_matrix")
#' weight_table <- inferCSN(example_matrix, verbose = TRUE)
#' weight_table_new <- crossweight(
#'   weight_table,
#'   matrix = t(example_matrix)
#' )
#' p1 <- network.heatmap(weight_table)
#' p2 <- network.heatmap(weight_table_new[, 1:3])
#' p1 + p2
#' }
crossweight <- function(
    weight_table,
    matrix,
    meta_data = NULL,
    lag = floor(ncol(matrix) / 5),
    min = ceiling(ncol(matrix) / 50),
    max = floor(ncol(matrix) / 12),
    symmetric_filter = FALSE,
    filter_thresh = 0) {
  if (!is.null(meta_data)) {
    matrix <- matrix[, rownames(meta_data$cell)]
  }
  weight_table$regulator <- as.character(weight_table$regulator)
  weight_table$target <- as.character(weight_table$target)
  weight_table$offset <- apply(weight_table, 1, cross_corr, matrix = matrix, lag = lag)

  weighted_score <- c()
  for (i in 1:nrow(weight_table)) {
    new <- score_offset(
      weight_table$weight[i],
      weight_table$offset[i],
      min = min,
      max = max,
      symmetric_filter = symmetric_filter
    )
    weighted_score <- c(weighted_score, new)
  }
  weight_table$weighted_score <- weighted_score

  # weight_table <- weight_table[abs(weight_table$weighted_score) > filter_thresh,]
  # weight_table <- weight_table[weight_table$weight > filter_thresh, ]
  weight_table <- weight_table[weight_table$offset > filter_thresh,]

  return(weight_table)
}

cross_corr <- function(
    grn_row,
    matrix,
    lag) {
  regulator <- grn_row[1]
  target <- grn_row[2]

  x <- stats::ccf(
    as.numeric(matrix[regulator, ]),
    as.numeric(matrix[target, ]),
    lag,
    plot = FALSE
  )

  df <- data.frame(lag = x$lag, cor = abs(x$acf))
  df <- df[order(df$cor, decreasing = TRUE), ]
  offset <- mean(df$lag[1:ceiling((2 / 3) * lag)])

  return(offset)
}

score_offset <- function(
    weight,
    offset,
    min = 2,
    max = 20,
    symmetric_filter = FALSE) {
  if (symmetric_filter) {
    offset <- abs(offset)
  }

  # if (offset <= min) {
  #   weight <- weight
  # } else if (offset >= max) {
  #   weight <- 0
  # } else {
  #   # Linear weighting scheme according to y = (-x / (max - min)) + 1
  #   score <- (-offset / (max - min)) + 1
  #   weight <- score * weight
  # }

  if (offset < 0 && weight < 0) {
    weight <- weight
  } else if (offset > 0 && weight > 0) {
    weight <- weight
  } else {
    weight <- 0
  }

  return(weight)
}

#' estimates min and max values for crossweighting for now assumes uniform cell density across
#' pseudotime/only considers early time this needs to be refined if it's to be useful...
#'
#' @param matrix matrix
#' @param meta_data meta_data
#' @param pseudotime_min pseudotime_min
#' @param pseudotime_max pseudotime_max
#'
#' @return Params list
#'
#' @export
crossweight_params <- function(
    matrix,
    meta_data,
    pseudotime_min = 0.005,
    pseudotime_max = 0.01) {
  matrix <- matrix[, rownames(meta_data$cell)]
  min <- nrow(meta_data$cell[meta_data$cell$pseudotime < pseudotime_min, ])
  max <- nrow(meta_data$cell[meta_data$cell$pseudotime < pseudotime_max, ])

  params <- list(
    matrix = matrix,
    min = min,
    max = max
  )

  return(params)
}
