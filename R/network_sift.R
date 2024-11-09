#' @title Sifting network
#'
#' @inheritParams inferCSN
#' @inheritParams network_format
#' @param matrix The expression matrix.
#' @param meta_data The meta data for cells or samples.
#' @param pseudotime_column The column of pseudotime.
#' @param method The method used for filter edges.
#' Could be choose \code{entropy} or \code{max}.
#' @param entropy_method If setting \code{method} to \code{entropy},
#'  could be choose \code{Shannon} or \code{Renyi} to compute entropy.
#' @param effective_entropy Default is \code{FALSE}.
#' Logical value, using effective entropy to filter weights or not.
#' @param shuffles Default is \code{100}.
#' The number of shuffles used to calculate the effective transfer entropy.
#' @param entropy_nboot Default is \code{300}.
#' The number of bootstrap replications for each direction of the estimated transfer entropy.
#' @param lag_value Default is \code{1}.
#' Markov order of gene expression values,
#' i.e. the number of lagged values affecting the current value of gene expression values.
#' @param entropy_p_value P value used to filter edges by entropy.
#'
#' @return Sifted network table
#' @export
#'
#' @examples
#' \dontrun{
#' data("example_matrix")
#' data("example_meta_data")
#' data("example_ground_truth")
#'
#' network_table <- inferCSN(example_matrix)
#' network_table_sifted <- network_sift(network_table)
#' network_table_sifted_entropy <- network_sift(
#'   network_table,
#'   matrix = example_matrix,
#'   meta_data = example_meta_data,
#'   pseudotime_column = "pseudotime",
#'   lag_value = 2,
#'   shuffles = 0,
#'   entropy_nboot = 0
#' )
#'
#' plot_network_heatmap(
#'   example_ground_truth[, 1:3],
#'   heatmap_title = "Ground truth",
#'   show_names = TRUE,
#'   rect_color = "gray70"
#' )
#' plot_network_heatmap(
#'   network_table,
#'   heatmap_title = "Raw",
#'   show_names = TRUE,
#'   rect_color = "gray70"
#' )
#' plot_network_heatmap(
#'   network_table_sifted,
#'   heatmap_title = "Filtered",
#'   show_names = TRUE,
#'   rect_color = "gray70"
#' )
#' plot_network_heatmap(
#'   network_table_sifted_entropy,
#'   heatmap_title = "Filtered by entropy",
#'   show_names = TRUE,
#'   rect_color = "gray70"
#' )
#'
#' calculate_metrics(
#'   network_table,
#'   example_ground_truth,
#'   plot = TRUE
#' )
#' calculate_metrics(
#'   network_table_sifted,
#'   example_ground_truth,
#'   plot = TRUE
#' )
#' calculate_metrics(
#'   network_table_sifted_entropy,
#'   example_ground_truth,
#'   plot = TRUE
#' )
#' }
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
    lag_value = 1,
    entropy_p_value = 0.05,
    cores = 1,
    verbose = TRUE) {
  method <- match.arg(method)
  if (method != "max") {
    entropy_method <- match.arg(entropy_method)
    if (is.null(matrix) | is.null(meta_data) | is.null(pseudotime_column)) {
      log_message(
        "Parameters: 'matrix', 'meta_data' and 'pseudotime_column' not all provide, setting 'method' to 'max'.",
        message_type = "warning",
        verbose = verbose
      )
      method <- "max"
      return(weight_sift(network_table))
    }
    if (!(pseudotime_column %in% colnames(meta_data))) {
      log_message(
        "Parameters: 'pseudotime_column' not in meta data provided, setting 'method' to 'max'.",
        message_type = "warning",
        verbose = verbose
      )
      method <- "max"
      return(weight_sift(network_table))
    }

    samples <- intersect(rownames(meta_data), rownames(matrix))
    if (is.null(samples)) {
      log_message(
        "No intersect samples in matrix and meta data, setting 'method' to 'max'.",
        message_type = "warning",
        verbose = verbose
      )
      method <- "max"
      return(weight_sift(network_table))
    }
  }

  if (method == "max") {
    return(weight_sift(network_table))
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
      log_message(
        "Parameter: 'effective_entropy == FALSE' and 'shuffles != 0', setting 'shuffles == 0'.",
        verbose = verbose
      )
      shuffles <- 0
    }
  } else {
    if (shuffles <= 10) {
      log_message(
        "Parameter: 'shuffles' is too small, setting 'shuffles == 10'.",
        verbose = verbose
      )
      shuffles <- 10
    }
  }

  transfer_entropy_table <- parallelize_fun(
    unique_pairs,
    cores = cores,
    verbose = verbose,
    fun = function(x) {
      result <- suppressWarnings(
        RTransferEntropy::transfer_entropy(
          matrix[, x[1]],
          matrix[, x[2]],
          lx = lag_value,
          ly = lag_value,
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
  ) |>
    purrr::list_rbind()

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
  ) |> weight_sift()

  network_table <- merge(
    network_table,
    transfer_entropy_table,
    by = c("regulator", "target")
  )
  network_table <- network_table[, c("regulator", "target", "weight")]

  return(network_table)
}
