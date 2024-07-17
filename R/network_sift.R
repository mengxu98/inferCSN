.weight_sift <- function(table) {
  table <- table[, 1:3]
  raw_rownames <- colnames(table)
  colnames(table) <- c("x", "y", "v")
  table$edge <- paste(
    table$x,
    table$y,
    sep = "_"
  )
  rownames(table) <- table$edge

  table_new <- data.frame(
    edge = paste(
      table$y,
      table$x,
      sep = "_"
    ),
    v = table$v
  )
  rownames(table_new) <- table_new$edge

  common_edges <- intersect(rownames(table), rownames(table_new))
  if (length(common_edges) == 0) {
    table <- table[, 1:3]
    colnames(table) <- raw_rownames
    rownames(table) <- NULL
    return(table)
  }
  table_common_edges <- table[common_edges, ]
  table_new <- table_new[common_edges, ]
  table_common_edges$v_new <- table_new$v
  table_common_edges <- dplyr::filter(table_common_edges, abs(v) < abs(v_new))
  table <- table[setdiff(rownames(table), rownames(table_common_edges)), 1:3]
  colnames(table) <- raw_rownames
  rownames(table) <- NULL

  return(table)
}

#' @title network_sift
#'
#' @inheritParams inferCSN
#' @param network_table network_table
#' @param matrix The expression matrix.
#' @param meta_data The meta data for cells or samples.
#' @param pseudotime_column The column of pseudotime.
#' @param method method The method used for filter edges. Could be choose \code{"entropy"} or \code{"max"}.
#' @param entropy_method If setting \code{'method'} to \code{'entropy'},
#'  could be choose \code{"Shannon"} or \code{"Renyi"} to compute entropy.
#' @param effective_entropy Logical value, using effective entropy to filter weights or not.
#' @param shuffles The number of shuffles used to calculate the effective transfer entropy. Default is \code{'shuffles'} = 100.
#' @param entropy_nboot entropy_nboot
#' @param history_length history_length
#' @param entropy_p_value P value.
#'
#' @return Filtered network table
#' @export
#'
#' @examples
#' data("example_matrix")
#' data("example_ground_truth")
#' network_table <- inferCSN(example_matrix)
#' network_table_filtered <- network_sift(network_table)
#' data("example_meta_data")
#' network_table_filtered_entropy <- network_sift(
#'   network_table,
#'   matrix = example_matrix,
#'   meta_data = example_meta_data,
#'   pseudotime_column = "pseudotime",
#'   history_length = 2,
#'   shuffles = 0,
#'   entropy_nboot = 0
#' )
#'
#' network.heatmap(
#'   example_ground_truth[, 1:3],
#'   heatmap_title = "Ground truth",
#'   show_names = TRUE,
#'   rect_color = "gray70"
#' )
#' network.heatmap(
#'   network_table,
#'   heatmap_title = "Raw",
#'   show_names = TRUE,
#'   rect_color = "gray70"
#' )
#' network.heatmap(
#'   network_table_filtered,
#'   heatmap_title = "Filtered",
#'   show_names = TRUE,
#'   rect_color = "gray70"
#' )
#' network.heatmap(
#'   network_table_filtered_entropy,
#'   heatmap_title = "Filtered by entropy",
#'   show_names = TRUE,
#'   rect_color = "gray70"
#' )
#'
#' calculate_auc(
#'   network_table,
#'   example_ground_truth,
#'   plot = TRUE
#' )
#' calculate_auc(
#'   network_table_filtered,
#'   example_ground_truth,
#'   plot = TRUE
#' )
#' calculate_auc(
#'   network_table_filtered_entropy,
#'   example_ground_truth,
#'   plot = TRUE
#' )
network_sift <- function(
    network_table,
    matrix = NULL,
    meta_data = NULL,
    pseudotime_column = NULL,
    method = c("entropy", "max"),
    entropy_method = c("Shannon", "Renyi"),
    effective_entropy = FALSE,
    shuffles = 100,
    entropy_nboot = 300,
    history_length = 1,
    entropy_p_value = 0.05,
    cores = 1,
    verbose = TRUE) {
  method <- match.arg(method)
  if (method != "max") {
    entropy_method <- match.arg(entropy_method)
    if (is.null(matrix) | is.null(meta_data) | is.null(pseudotime_column)) {
      if (verbose) {
        message(
          "Parameters: 'matrix', 'meta_data' and 'pseudotime_column' not all provide, setting 'method' to 'max'."
        )
      }
      method <- "max"
      return(.weight_sift(network_table))
    }
    if (!(pseudotime_column %in% colnames(meta_data))) {
      if (verbose) {
        message(
          "Parameters: 'pseudotime_column' not in meta data provided, setting 'method' to 'max'."
        )
      }
      method <- "max"
      return(.weight_sift(network_table))
    }

    samples <- intersect(rownames(meta_data), rownames(matrix))
    if (is.null(samples)) {
      if (verbose) {
        message(
          "No intersect samples in matrix and meta data, setting 'method' to 'max'."
        )
      }
      method <- "max"
      return(.weight_sift(network_table))
    }
  }

  if (method == "max") {
    return(.weight_sift(network_table))
  }

  meta_data <- meta_data[samples, ]
  meta_data <- meta_data[order(
    meta_data[, pseudotime_column],
    decreasing = FALSE
  ), ]

  genes <- unique(c(network_table$regulator, network_table$target))

  matrix <- matrix[rownames(meta_data), genes]

  pairs <- expand.grid(
    regulator = colnames(matrix),
    target = colnames(matrix)
  )
  pairs <- pairs[pairs$regulator != pairs$target, ]
  pairs <- as.data.frame(pairs)

  pairs <- purrr::map(
    seq_len(nrow(pairs)),
    .f = function(x) {
      as.character(c(pairs[x, 1], pairs[x, 2]))
    }
  )

  unique_pairs <- list()
  for (pair in pairs) {
    ordered_pair <- sort(pair)

    found <- FALSE
    for (up in unique_pairs) {
      if (all(up == ordered_pair)) {
        found <- TRUE
        break
      }
    }

    if (!found) {
      unique_pairs <- append(unique_pairs, list(ordered_pair))
    }
  }

  if (!effective_entropy) {
    if (shuffles != 0) {
      if (verbose) {
        message(
          "Parameter: 'effective_entropy == FALSE' and 'shuffles != 0', setting 'shuffles == 0'."
        )
      }
      shuffles <- 0
    }
  } else {
    if (shuffles <= 10) {
      if (verbose) {
        message(
          "Parameter: 'shuffles' is too small, setting 'shuffles == 10'."
        )
      }
      shuffles <- 10
    }
  }

  transfer_entropy_list <- parallelize_fun(
    unique_pairs,
    cores = cores,
    verbose = verbose,
    fun = function(x) {
      result <- suppressWarnings(
        RTransferEntropy::transfer_entropy(
          matrix[, x[1]],
          matrix[, x[2]],
          lx = history_length,
          ly = history_length,
          entropy = entropy_method,
          shuffles = shuffles,
          nboot = entropy_nboot,
          quiet = TRUE
        )
      )

      result <- coef(result)
      if (effective_entropy) {
        data.frame(
          "regulator" = x[1],
          "target" = x[2],
          "entropy" = result[1, 2],
          "entropy_contrary" = result[2, 2],
          "P_value" = result[1, 4],
          "P_value_contrary" = result[2, 4]
        )
      } else {
        data.frame(
          "regulator" = x[1],
          "target" = x[2],
          "entropy" = result[1, 1],
          "entropy_contrary" = result[2, 1],
          "P_value" = result[1, 4],
          "P_value_contrary" = result[2, 4]
        )
      }
    }
  )

  transfer_entropy_table <- purrr::map_dfr(
    transfer_entropy_list,
    .f = function(x) {
      x
    }
  )

  if (entropy_nboot > 1) {
    transfer_entropy_table <- dplyr::filter(
      transfer_entropy_table,
      P_value <= entropy_p_value & P_value_contrary <= entropy_p_value
    )
  }

  transfer_entropy_table_contrary <- transfer_entropy_table[, c(2, 1, 4)]
  colnames(transfer_entropy_table_contrary) <- c("regulator", "target", "entropy")
  transfer_entropy_table <- rbind(
    transfer_entropy_table[, c(1, 2, 3)],
    transfer_entropy_table_contrary
  )

  transfer_entropy_table <- .weight_sift(transfer_entropy_table)

  network_table <- merge(
    network_table,
    transfer_entropy_table,
    by = c("regulator", "target")
  )
  network_table <- network_table[, c("regulator", "target", "weight")]

  return(network_table)
}
