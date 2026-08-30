#' @title Calculate network metrics
#'
#' @md
#' @param network_table Predicted edge table.
#' @param ground_truth Ground-truth edge table.
#' @param metric_type Metrics to compute.
#' @param return_plot Whether to return a plot.
#' @param line_color,line_width Plot line style.
#' @param tf_edges Whether to restrict regulators to transcription factors.
#' @return A list with metrics and an optional plot.
#' @export
calculate_metrics <- function(
  network_table,
  ground_truth,
  metric_type = c(
    "all",
    "all_no_epr",
    "auc",
    "auroc",
    "auprc",
    "epr",
    "precision",
    "recall",
    "f1",
    "accuracy",
    "si",
    "ji"
  ),
  return_plot = FALSE,
  tf_edges = FALSE,
  line_color = "#1563cc",
  line_width = 1
) {
  metric_type <- match.arg(metric_type)

  result <- switch(
    EXPR = metric_type,
    "all" = {
      pred_data <- prepare_binary_predictions(
        network_table,
        ground_truth,
        tf_edges = tf_edges
      )
      calculate_all_metrics_from_pred_data(
        pred_data,
        network_table,
        ground_truth,
        include_epr = TRUE,
        tf_edges = tf_edges,
        return_plot = return_plot,
        line_color = line_color,
        line_width = line_width
      )
    },
    "all_no_epr" = {
      pred_data <- prepare_binary_predictions(
        network_table,
        ground_truth,
        tf_edges = tf_edges
      )
      calculate_all_metrics_from_pred_data(
        pred_data,
        network_table,
        ground_truth,
        include_epr = FALSE,
        tf_edges = tf_edges,
        return_plot = return_plot,
        line_color = line_color,
        line_width = line_width
      )
    },
    "auc" = calculate_auc(
      network_table = network_table,
      ground_truth = ground_truth,
      return_plot = return_plot,
      line_color = line_color,
      line_width = line_width,
      tf_edges = tf_edges
    ),
    "auroc" = calculate_auroc(
      network_table = network_table,
      ground_truth = ground_truth,
      return_plot = return_plot,
      line_color = line_color,
      line_width = line_width,
      tf_edges = tf_edges
    ),
    "auprc" = calculate_auprc(
      network_table = network_table,
      ground_truth = ground_truth,
      return_plot = return_plot,
      line_color = line_color,
      line_width = line_width,
      tf_edges = tf_edges
    ),
    "epr" = calculate_epr(network_table, ground_truth, tf_edges = tf_edges),
    "precision" = calculate_precision(
      network_table,
      ground_truth,
      tf_edges = tf_edges
    ),
    "recall" = calculate_recall(
      network_table,
      ground_truth,
      tf_edges = tf_edges
    ),
    "f1" = calculate_f1(
      network_table,
      ground_truth,
      tf_edges = tf_edges
    ),
    "accuracy" = calculate_accuracy(
      network_table,
      ground_truth,
      tf_edges = tf_edges
    ),
    "si" = calculate_si(
      network_table,
      ground_truth,
      tf_edges = tf_edges
    ),
    "ji" = calculate_ji(
      network_table,
      ground_truth,
      tf_edges = tf_edges
    )
  )

  result
}

safe_metric_divide <- function(numerator, denominator, default = 0) {
  if (is.na(denominator) || denominator == 0) {
    return(default)
  }
  numerator / denominator
}


calculate_auroc_from_gold <- function(
  gold,
  return_plot = FALSE,
  line_color = "#1563cc",
  line_width = 1
) {
  auroc <- as.numeric(sprintf(
    "%0.3f",
    compute_rank_auroc(gold$weight, gold$label)
  ))
  result <- data.frame(Metric = "AUROC", Value = auroc)

  if (return_plot) {
    auc_curves <- precrec::evalmod(scores = gold$weight, labels = gold$label)
    return(list(
      metrics = result,
      plot = plot_auroc(auc_curves, auroc, line_color, line_width)
    ))
  }

  list(metrics = result)
}

calculate_auprc_from_gold <- function(
  gold,
  return_plot = FALSE,
  line_color = "#1563cc",
  line_width = 1
) {
  auprc <- as.numeric(sprintf(
    "%0.3f",
    compute_rank_auprc(gold$weight, gold$label)
  ))
  result <- data.frame(Metric = "AUPRC", Value = auprc)

  if (return_plot) {
    auc_curves <- precrec::evalmod(scores = gold$weight, labels = gold$label)
    return(list(
      metrics = result,
      plot = plot_auprc(auc_curves, auprc, line_color, line_width)
    ))
  }

  list(metrics = result)
}

calculate_auc_from_pred_data <- function(
  pred_data,
  return_plot = FALSE,
  line_color = "#1563cc",
  line_width = 1,
  tag_levels = "A"
) {
  gold <- pred_data$gold
  auroc <- as.numeric(sprintf(
    "%0.3f",
    compute_rank_auroc(gold$weight, gold$label)
  ))
  auprc <- as.numeric(sprintf(
    "%0.3f",
    compute_rank_auprc(gold$weight, gold$label)
  ))
  metrics <- data.frame(
    Metric = c("AUROC", "AUPRC"),
    Value = c(auroc, auprc)
  )

  if (return_plot) {
    auc_curves <- precrec::evalmod(scores = gold$weight, labels = gold$label)
    final_plot <- patchwork::wrap_plots(
      list(
        plot_auroc(auc_curves, auroc, line_color, line_width),
        plot_auprc(auc_curves, auprc, line_color, line_width)
      ),
      ncol = 2
    ) +
      patchwork::plot_annotation(tag_levels = tag_levels)
    return(list(metrics = metrics, plot = final_plot))
  }

  list(metrics = metrics)
}

binary_counts_from_pred_data <- function(pred_data) {
  counts <- as.vector(
    table(
      pred_data$predictor_binary,
      pred_data$true_label,
      dnn = c("Predicted", "Actual")
    )
  )
  names(counts) <- c("TN", "FN", "FP", "TP")
  counts
}

calculate_all_metrics_from_pred_data <- function(
  pred_data,
  network_table,
  ground_truth,
  include_epr = TRUE,
  tf_edges = FALSE,
  return_plot = FALSE,
  line_color = "#1563cc",
  line_width = 1
) {
  auc_result <- calculate_auc_from_pred_data(
    pred_data,
    return_plot = FALSE,
    line_color = line_color,
    line_width = line_width
  )
  epr_result <- if (isTRUE(include_epr)) {
    calculate_epr(network_table, ground_truth, tf_edges = tf_edges)
  } else {
    NULL
  }
  counts <- binary_counts_from_pred_data(pred_data)

  precision_raw <- safe_metric_divide(
    counts[["TP"]],
    counts[["TP"]] + counts[["FP"]],
    0
  )
  recall_raw <- safe_metric_divide(
    counts[["TP"]],
    counts[["TP"]] + counts[["FN"]],
    0
  )
  precision <- as.numeric(sprintf("%0.3f", precision_raw))
  recall <- as.numeric(sprintf("%0.3f", recall_raw))
  f1_raw <- if ((precision_raw + recall_raw) == 0) {
    0
  } else {
    2 * precision_raw * recall_raw / (precision_raw + recall_raw)
  }
  f1 <- as.numeric(sprintf("%0.3f", f1_raw))
  accuracy <- as.numeric(sprintf(
    "%0.3f",
    safe_metric_divide(counts[["TP"]] + counts[["TN"]], sum(counts), 0)
  ))

  overlap <- length(intersect(pred_data$pred_edge_ids, pred_data$true_edge_ids))
  union_n <- length(unique(c(pred_data$pred_edge_ids, pred_data$true_edge_ids)))
  ji <- if (union_n == 0) 0 else overlap / union_n
  ji <- as.numeric(sprintf("%0.3f", ji))

  metrics_df <- rbind(
    auc_result$metrics,
    if (!is.null(epr_result)) epr_result$metrics,
    data.frame(Metric = "Precision", Value = precision),
    data.frame(Metric = "Recall", Value = recall),
    data.frame(Metric = "F1", Value = f1),
    data.frame(Metric = "Accuracy", Value = accuracy),
    data.frame(Metric = "JI", Value = ji),
    data.frame(Metric = "SI", Value = overlap)
  )
  rownames(metrics_df) <- NULL

  if (return_plot) {
    plot_metrics <- metrics_df[metrics_df$Metric != "SI", , drop = FALSE]
    return(list(
      metrics = metrics_df,
      plot = plot_all_metrics(plot_metrics, line_color)
    ))
  }

  list(metrics = metrics_df)
}

normalize_ground_truth_edges <- function(ground_truth) {
  truth <- as.data.frame(ground_truth, stringsAsFactors = FALSE)
  if (ncol(truth) > 2) {
    truth <- truth[, 1:2, drop = FALSE]
  }
  names(truth)[1:2] <- c("regulator", "target")
  truth <- truth[!is.na(truth$regulator) & !is.na(truth$target), , drop = FALSE]
  truth <- truth[nzchar(truth$regulator) & nzchar(truth$target), , drop = FALSE]
  truth <- truth[truth$regulator != truth$target, , drop = FALSE]
  unique(truth)
}

normalize_signed_predicted_edges <- function(network_table) {
  pred <- as.data.frame(network_table, stringsAsFactors = FALSE)
  keep <- intersect(c("regulator", "target", "weight"), names(pred))
  if (length(keep) < 3) {
    stop("network_table must contain columns: regulator, target, weight")
  }
  pred <- pred[, keep, drop = FALSE]
  pred$weight <- suppressWarnings(as.numeric(pred$weight))
  pred <- pred[is.finite(pred$weight), , drop = FALSE]
  pred <- pred[!is.na(pred$regulator) & !is.na(pred$target), , drop = FALSE]
  pred <- pred[nzchar(pred$regulator) & nzchar(pred$target), , drop = FALSE]
  pred <- pred[pred$regulator != pred$target, , drop = FALSE]
  pred$abs_weight <- abs(pred$weight)
  pred <- pred[order(pred$abs_weight, decreasing = TRUE), , drop = FALSE]
  pred <- pred[!duplicated(pred[, c("regulator", "target")]), , drop = FALSE]
  rownames(pred) <- NULL
  pred
}

normalize_signed_ground_truth_edges <- function(ground_truth) {
  truth <- as.data.frame(ground_truth, stringsAsFactors = FALSE)
  if (ncol(truth) < 3) {
    stop(
      "ground_truth must contain a third column with edge sign/type information"
    )
  }
  truth <- truth[, 1:3, drop = FALSE]
  names(truth)[1:3] <- c("regulator", "target", "type")
  truth <- truth[
    !is.na(truth$regulator) & !is.na(truth$target) & !is.na(truth$type), ,
    drop = FALSE
  ]
  truth <- truth[nzchar(truth$regulator) & nzchar(truth$target), , drop = FALSE]
  truth <- truth[truth$regulator != truth$target, , drop = FALSE]
  truth$type <- as.character(truth$type)
  truth <- truth[truth$type %in% c("+", "-"), , drop = FALSE]
  unique(truth)
}

select_top_ranked_edges <- function(
  network_table,
  top_k,
  drop_zero = FALSE,
  restrict_edges = NULL
) {
  pred <- normalize_signed_predicted_edges(network_table)

  if (!is.null(restrict_edges)) {
    edge_ids <- paste(pred$regulator, pred$target, sep = "\r")
    pred <- pred[edge_ids %in% restrict_edges, , drop = FALSE]
  }

  if (drop_zero) {
    pred <- pred[pred$abs_weight > 0, , drop = FALSE]
  }

  if (!nrow(pred)) {
    return(pred)
  }

  maxk <- min(nrow(pred), top_k)
  if (maxk <= 0) {
    return(pred[0, , drop = FALSE])
  }
  pred[seq_len(maxk), , drop = FALSE]
}

select_top_ranked_from_normalized <- function(pred, top_k, drop_zero = FALSE) {
  if (!nrow(pred) || top_k <= 0) {
    return(pred[0, , drop = FALSE])
  }

  if (drop_zero) {
    pred <- pred[pred$abs_weight > 0, , drop = FALSE]
  }
  if (!nrow(pred)) {
    return(pred)
  }

  maxk <- min(nrow(pred), top_k)
  pred[seq_len(maxk), , drop = FALSE]
}

pairwise_summary <- function(values, metric_name, mad_name) {
  values <- values[is.finite(values)]
  data.frame(
    Metric = c(metric_name, mad_name),
    Value = c(
      if (length(values)) {
        as.numeric(sprintf("%0.3f", stats::median(values)))
      } else {
        NA_real_
      },
      if (length(values)) {
        as.numeric(sprintf("%0.3f", mean(abs(values - mean(values)))))
      } else {
        NA_real_
      }
    )
  )
}

format_metric_values <- function(values, digits = 3) {
  vapply(
    values,
    function(x) {
      if (is.na(x)) {
        NA_real_
      } else {
        as.numeric(sprintf(paste0("%0.", digits, "f"), x))
      }
    },
    numeric(1)
  )
}

compute_network_scores <- function(network_table, ground_truth, tf_edges = FALSE) {
  pred <- normalize_signed_predicted_edges(network_table)
  pred$weight <- abs(pred$weight)
  truth <- normalize_ground_truth_edges(ground_truth)

  gt_genes <- sort(unique(c(truth$regulator, truth$target)))
  if (!length(gt_genes)) {
    empty <- data.frame(
      regulator = character(0),
      target = character(0),
      weight = numeric(0),
      label = integer(0),
      stringsAsFactors = FALSE
    )
    return(list(
      gold = empty,
      pred_edge_ids = character(0),
      true_edge_ids = character(0),
      predictor_binary = factor(integer(0), levels = c(0, 1)),
      true_label = factor(integer(0), levels = c(0, 1))
    ))
  }

  gold <- if (exists("prepare_calculate_metrics", mode = "function") && !isTRUE(tf_edges)) {
    prepare_calculate_metrics(pred, truth)
  } else {
    if (isTRUE(tf_edges)) {
      tf_genes <- sort(unique(truth$regulator))
      universe <- expand.grid(
        regulator = tf_genes,
        target = gt_genes,
        KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
      )
    } else {
      universe <- expand.grid(
        regulator = gt_genes,
        target = gt_genes,
        KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
      )
    }
    universe <- universe[universe$regulator != universe$target, , drop = FALSE]

    truth_edge_ids <- paste(
      truth$regulator,
      truth$target,
      sep = "\n"
    )
    pred_edge_ids <- paste(
      pred$regulator,
      pred$target,
      sep = "\n"
    )
    pred_lookup <- stats::setNames(pred$weight, pred_edge_ids)

    universe_edge_ids <- paste(
      universe$regulator,
      universe$target,
      sep = "\n"
    )
    universe$weight <- pred_lookup[universe_edge_ids]
    universe$weight[is.na(universe$weight)] <- 0
    universe$label <- as.integer(universe_edge_ids %in% truth_edge_ids)
    universe
  }

  threshold_result <- select_best_binary_threshold(
    gold$weight,
    gold$label
  )
  predictor_binary <- factor(
    ifelse(threshold_result$predicted_positive, 1, 0),
    levels = c(0, 1)
  )
  true_label <- factor(gold$label, levels = c(0, 1))

  pred_binary_edges <- gold[threshold_result$predicted_positive, c("regulator", "target"), drop = FALSE]
  pred_edge_ids <- paste(pred_binary_edges$regulator, pred_binary_edges$target, sep = "-")
  true_edge_ids <- paste(truth$regulator, truth$target, sep = "-")

  list(
    gold = gold,
    pred_edge_ids = pred_edge_ids,
    true_edge_ids = true_edge_ids,
    predictor_binary = predictor_binary,
    true_label = true_label
  )
}

compute_rank_auroc <- function(scores, labels) {
  positive <- labels == 1
  negative <- labels == 0
  n_pos <- as.numeric(sum(positive))
  n_neg <- as.numeric(sum(negative))

  if (n_pos == 0 || n_neg == 0) {
    return(NA_real_)
  }

  ranks <- rank(scores, ties.method = "average")
  safe_metric_divide(
    sum(ranks[positive]) - (n_pos * (n_pos + 1) / 2),
    n_pos * n_neg,
    default = NA_real_
  )
}

precision_recall_curve <- function(scores, labels) {
  ord <- order(scores, decreasing = TRUE)
  scores_ord <- scores[ord]
  labels_ord <- labels[ord]

  tp <- cumsum(labels_ord == 1)
  fp <- cumsum(labels_ord == 0)
  threshold_idx <- c(which(diff(scores_ord) != 0), length(scores_ord))
  tps <- tp[threshold_idx]
  fps <- fp[threshold_idx]
  precision <- ifelse((tps + fps) > 0, tps / (tps + fps), 0)

  if (length(tps) == 0 || tail(tps, 1) == 0) {
    recall <- rep(1, length(tps))
  } else {
    recall <- tps / tail(tps, 1)
  }

  rev_idx <- rev(seq_along(precision))
  list(
    precision = c(precision[rev_idx], 1),
    recall = c(recall[rev_idx], 0)
  )
}

compute_rank_auprc <- function(scores, labels) {
  n_pos <- sum(labels == 1)
  if (n_pos == 0) {
    return(NA_real_)
  }

  curve <- precision_recall_curve(scores, labels)
  abs(sum(
    diff(curve$recall) *
      (head(curve$precision, -1) + tail(curve$precision, -1)) /
      2
  ))
}

select_best_binary_threshold <- function(scores, labels) {
  thresholds <- sort(unique(scores), decreasing = TRUE)
  if (length(thresholds) == 0) {
    thresholds <- 0
  }

  best <- list(
    threshold = thresholds[[1]],
    youden = -Inf,
    predicted_positive = rep(FALSE, length(scores))
  )

  for (threshold in thresholds) {
    predicted_positive <- scores >= threshold
    tp <- sum(predicted_positive & labels == 1)
    fp <- sum(predicted_positive & labels == 0)
    tn <- sum(!predicted_positive & labels == 0)
    fn <- sum(!predicted_positive & labels == 1)

    sensitivity <- safe_metric_divide(tp, tp + fn, default = NA_real_)
    specificity <- safe_metric_divide(tn, tn + fp, default = NA_real_)
    youden <- sensitivity + specificity - 1
    if (is.na(youden)) {
      youden <- -Inf
    }

    if (youden > best$youden) {
      best <- list(
        threshold = threshold,
        youden = youden,
        predicted_positive = predicted_positive
      )
    }
  }

  best
}

prepare_binary_predictions <- function(network_table, ground_truth, tf_edges = FALSE) {
  compute_network_scores(network_table, ground_truth, tf_edges = tf_edges)
}

plot_auroc <- function(auc_curves, aurco, line_color, line_width) {
  subset(fortify(auc_curves), curvetype == "ROC") |>
    ggplot(aes(x = x, y = y)) +
    geom_line(color = line_color, linewidth = line_width) +
    geom_abline(
      slope = 1,
      intercept = 0,
      color = line_color,
      linetype = "dotted",
      linewidth = line_width
    ) +
    labs(
      title = paste("AUROC:", aurco),
      x = "False positive rate",
      y = "True positive rate"
    ) +
    theme_bw() +
    coord_fixed()
}

plot_auprc <- function(auc_curves, auprc, line_color, line_width) {
  subset(fortify(auc_curves), curvetype == "PRC") |>
    ggplot(aes(x = x, y = y)) +
    geom_line(color = line_color, linewidth = line_width) +
    labs(
      title = paste("AUPRC:", auprc),
      x = "Recall",
      y = "Precision"
    ) +
    theme_bw() +
    coord_fixed()
}

#' @title Calculate AUROC Metric
#' @description Calculates AUROC metric with optional visualization
#' @param network_table A data frame of predicted network structure
#' @param ground_truth A data frame of ground truth network
#' @param return_plot Logical value indicating whether to generate plot
#' @param line_color Color for plot lines
#' @param tf_edges Whether to restrict candidate edges to TF-to-gene.
#'   Default is `FALSE`.
#' @param line_width Width for plot lines
#' @return A list containing metric and optional plot
#' @export
calculate_auroc <- function(
  network_table,
  ground_truth,
  return_plot = FALSE,
  line_color = "#1563cc",
  line_width = 1,
  tf_edges = FALSE
) {
  pred_data <- prepare_binary_predictions(
    network_table,
    ground_truth,
    tf_edges = tf_edges
  )
  calculate_auroc_from_gold(
    pred_data$gold,
    return_plot = return_plot,
    line_color = line_color,
    line_width = line_width
  )
}

#' @title Calculate AUPRC Metric
#' @description Calculates AUPRC metric with optional visualization
#' @param network_table A data frame of predicted network structure
#' @param ground_truth A data frame of ground truth network
#' @param return_plot Logical value indicating whether to generate plot
#' @param line_color Color for plot lines
#' @param tf_edges Whether to restrict candidate edges to TF-to-gene.
#'   Default is `FALSE`.
#' @param line_width Width for plot lines
#' @return A list containing metric and optional plot
#' @export
calculate_auprc <- function(
  network_table,
  ground_truth,
  return_plot = FALSE,
  line_color = "#1563cc",
  line_width = 1,
  tf_edges = FALSE
) {
  pred_data <- prepare_binary_predictions(
    network_table,
    ground_truth,
    tf_edges = tf_edges
  )
  calculate_auprc_from_gold(
    pred_data$gold,
    return_plot = return_plot,
    line_color = line_color,
    line_width = line_width
  )
}

#' @title Calculate Early Precision Ratio
#' @description Calculates the Early Precision Ratio (EPR) on the fixed candidate universe.
#' EPR compares the precision among the top-ranked predicted edges to the
#' precision expected from a random predictor over the same candidate edge pool.
#' @param network_table A data frame of predicted network structure
#' @param ground_truth A data frame of ground truth network
#' @param tf_edges Whether to restrict the candidate universe to regulator-to-gene edges
#' @return A list containing the metric
#' @export
calculate_epr <- function(
  network_table,
  ground_truth,
  tf_edges = FALSE
) {
  truth <- normalize_ground_truth_edges(ground_truth)
  if (!nrow(truth)) {
    value <- NA_real_
    value <- format_metric_values(value)
    return(list(metrics = data.frame(Metric = "EPR", Value = value)))
  }

  universe_genes <- sort(unique(c(truth$regulator, truth$target)))
  n_genes <- length(universe_genes)
  true_edges <- unique(paste(truth$regulator, truth$target, sep = "\r"))
  n_true <- length(true_edges)
  n_possible <- if (isTRUE(tf_edges)) {
    length(unique(truth$regulator)) * max(n_genes - 1L, 0L)
  } else {
    n_genes * max(n_genes - 1L, 0L)
  }

  if (n_true == 0 || n_possible == 0) {
    value <- NA_real_
  } else {
    pred <- normalize_signed_predicted_edges(network_table)
    if (tf_edges) {
      tf_genes <- unique(truth$regulator)
      pred <- pred[
        pred$regulator %in% tf_genes &
          pred$target %in% universe_genes, ,
        drop = FALSE
      ]
      pred <- select_top_ranked_from_normalized(pred, top_k = n_true)
    } else {
      pred <- select_top_ranked_from_normalized(pred, top_k = n_true)
    }

    if (!nrow(pred)) {
      value <- 0
    } else {
      pred_ids <- paste(pred$regulator, pred$target, sep = "\r")
      eprec <- length(intersect(pred_ids, true_edges)) / n_true
      value <- safe_metric_divide(
        eprec,
        n_true / n_possible,
        default = NA_real_
      )
    }
  }

  value <- format_metric_values(value)
  list(metrics = data.frame(Metric = "EPR", Value = value))
}

#' @title Calculate AUC Metrics
#' @description Calculates AUROC and AUPRC metrics with optional visualization
#' @param network_table A data frame of predicted network structure
#' @param ground_truth A data frame of ground truth network
#' @param return_plot Logical value indicating whether to generate plots
#' @param line_color Color for plot lines
#' @param line_width Width for plot lines
#' @param tag_levels Tag levels for plot annotations
#' @param tf_edges Whether to restrict candidate edges to TF-to-gene.
#'   Default is `FALSE`.
#' @return A list containing metrics and optional plots
#' @export
calculate_auc <- function(
  network_table,
  ground_truth,
  return_plot = FALSE,
  line_color = "#1563cc",
  line_width = 1,
  tag_levels = "A",
  tf_edges = FALSE
) {
  pred_data <- prepare_binary_predictions(
    network_table,
    ground_truth,
    tf_edges = tf_edges
  )
  calculate_auc_from_pred_data(
    pred_data,
    return_plot = return_plot,
    line_color = line_color,
    line_width = line_width,
    tag_levels = tag_levels
  )
}

metric_from_binary <- function(network_table, ground_truth, tf_edges = FALSE) {
  pred_data <- prepare_binary_predictions(
    network_table,
    ground_truth,
    tf_edges = tf_edges
  )
  binary_counts_from_pred_data(pred_data)
}

#' @title Calculate Precision
#' @description Calculates precision (positive predictive value) for predicted edges
#'   against a ground truth network after optimal thresholding.
#' @param network_table A data frame of predicted network structure.
#' @param ground_truth A data frame of ground truth network.
#' @param tf_edges Whether to restrict the candidate edge universe to
#'   TF-to-gene edges. Default is `FALSE`.
#' @return A list with a `metrics` data frame containing Precision.
#' @export
calculate_precision <- function(network_table, ground_truth, tf_edges = FALSE) {
  counts <- metric_from_binary(network_table, ground_truth, tf_edges = tf_edges)
  precision <- as.numeric(sprintf(
    "%0.3f",
    safe_metric_divide(counts[["TP"]], counts[["TP"]] + counts[["FP"]], 0)
  ))
  list(metrics = data.frame(Metric = "Precision", Value = precision))
}

#' @title Calculate Recall
#' @description Calculates recall (true positive rate) for predicted edges
#'   against a ground truth network after optimal thresholding.
#' @param network_table A data frame of predicted network structure.
#' @param ground_truth A data frame of ground truth network.
#' @param tf_edges Whether to restrict the candidate edge universe to
#'   TF-to-gene edges. Default is `FALSE`.
#' @return A list with a `metrics` data frame containing Recall.
#' @export
calculate_recall <- function(network_table, ground_truth, tf_edges = FALSE) {
  counts <- metric_from_binary(network_table, ground_truth, tf_edges = tf_edges)
  recall <- as.numeric(sprintf(
    "%0.3f",
    safe_metric_divide(counts[["TP"]], counts[["TP"]] + counts[["FN"]], 0)
  ))
  list(metrics = data.frame(Metric = "Recall", Value = recall))
}

#' @title Calculate F1 Score
#' @description Calculates the F1 score (harmonic mean of precision and recall)
#'   for predicted edges against a ground truth network after optimal thresholding.
#' @param network_table A data frame of predicted network structure.
#' @param ground_truth A data frame of ground truth network.
#' @param tf_edges Whether to restrict the candidate edge universe to
#'   TF-to-gene edges. Default is `FALSE`.
#' @return A list with a `metrics` data frame containing F1.
#' @export
calculate_f1 <- function(network_table, ground_truth, tf_edges = FALSE) {
  counts <- metric_from_binary(network_table, ground_truth, tf_edges = tf_edges)
  precision <- safe_metric_divide(
    counts[["TP"]],
    counts[["TP"]] + counts[["FP"]],
    0
  )
  recall <- safe_metric_divide(
    counts[["TP"]],
    counts[["TP"]] + counts[["FN"]],
    0
  )
  f1 <- if ((precision + recall) == 0) {
    0
  } else {
    2 * precision * recall / (precision + recall)
  }
  f1 <- as.numeric(sprintf("%0.3f", f1))
  list(metrics = data.frame(Metric = "F1", Value = f1))
}

#' @title Calculate Accuracy
#' @description Calculates overall classification accuracy for predicted edges
#'   against a ground truth network after optimal thresholding.
#' @param network_table A data frame of predicted network structure.
#' @param ground_truth A data frame of ground truth network.
#' @param tf_edges Whether to restrict the candidate edge universe to
#'   TF-to-gene edges. Default is `FALSE`.
#' @return A list with a `metrics` data frame containing Accuracy.
#' @export
calculate_accuracy <- function(network_table, ground_truth, tf_edges = FALSE) {
  counts <- metric_from_binary(network_table, ground_truth, tf_edges = tf_edges)
  accuracy <- as.numeric(sprintf(
    "%0.3f",
    safe_metric_divide(counts[["TP"]] + counts[["TN"]], sum(counts), 0)
  ))
  list(metrics = data.frame(Metric = "Accuracy", Value = accuracy))
}

#' @title Calculate Set Intersection
#' @description Calculates the Set Intersection (SI) — the number of predicted
#'   edges that exactly match true edges in the ground truth network.
#' @param network_table A data frame of predicted network structure.
#' @param ground_truth A data frame of ground truth network.
#' @param tf_edges Whether to restrict the candidate edge universe to
#'   TF-to-gene edges. Default is `FALSE`.
#' @return A list with a `metrics` data frame containing SI.
#' @export
calculate_si <- function(network_table, ground_truth, tf_edges = FALSE) {
  pred_data <- prepare_binary_predictions(
    network_table,
    ground_truth,
    tf_edges = tf_edges
  )
  si <- length(intersect(pred_data$pred_edge_ids, pred_data$true_edge_ids))
  list(metrics = data.frame(Metric = "SI", Value = si))
}

#' @title Calculate Jaccard Index
#' @description Calculates the Jaccard Index (JI) — the size of the intersection
#'   divided by the size of the union of predicted and true edge sets.
#' @param network_table A data frame of predicted network structure.
#' @param ground_truth A data frame of ground truth network.
#' @param tf_edges Whether to restrict the candidate edge universe to
#'   TF-to-gene edges. Default is `FALSE`.
#' @return A list with a `metrics` data frame containing JI.
#' @export
calculate_ji <- function(network_table, ground_truth, tf_edges = FALSE) {
  pred_data <- prepare_binary_predictions(
    network_table,
    ground_truth,
    tf_edges = tf_edges
  )
  overlap <- length(intersect(pred_data$pred_edge_ids, pred_data$true_edge_ids))
  union_n <- length(unique(c(pred_data$pred_edge_ids, pred_data$true_edge_ids)))
  ji <- if (union_n == 0) 0 else overlap / union_n
  ji <- as.numeric(sprintf("%0.3f", ji))
  list(metrics = data.frame(Metric = "JI", Value = ji))
}

#' @title Calculate Signed Early Precision Ratio
#' @description Calculates activation and inhibitory early precision ratios using a signed ground truth.
#' These metrics evaluate positive and negative regulatory edges separately,
#' using the same random-baseline normalization as EPR.
#' @param network_table A data frame of predicted network structure
#' @param ground_truth A data frame with columns regulator, target and type (`+` or `-`)
#' @return A list containing activation and inhibitory EPR metrics
#' @export
calculate_signed_epr <- function(network_table, ground_truth) {
  truth <- normalize_signed_ground_truth_edges(ground_truth)
  activation_ids <- paste(
    truth$regulator[truth$type == "+"],
    truth$target[truth$type == "+"],
    sep = "\r"
  )
  inhibitory_ids <- paste(
    truth$regulator[truth$type == "-"],
    truth$target[truth$type == "-"],
    sep = "\r"
  )
  all_genes <- unique(c(truth$regulator, truth$target))
  n_possible <- length(all_genes) * max(length(all_genes) - 1, 0)
  pred <- normalize_signed_predicted_edges(network_table)

  signed_epr_one <- function(true_ids, opposite_ids) {
    k <- length(unique(true_ids))
    if (k == 0 || n_possible == 0) {
      return(NA_real_)
    }

    pred_ids <- paste(pred$regulator, pred$target, sep = "\r")
    candidates <- pred[!(pred_ids %in% opposite_ids), , drop = FALSE]
    if (!nrow(candidates)) {
      return(NA_real_)
    }

    maxk <- min(nrow(candidates), k)
    edge_weight_topk <- candidates$abs_weight[[maxk]]
    non_zero_min <- if (any(candidates$abs_weight > 0)) {
      min(candidates$abs_weight[candidates$abs_weight > 0])
    } else {
      0
    }
    best_val <- max(non_zero_min, edge_weight_topk)
    selected <- candidates[candidates$abs_weight >= best_val, , drop = FALSE]
    if (!nrow(selected)) {
      return(NA_real_)
    }

    selected_ids <- paste(selected$regulator, selected$target, sep = "\r")
    early_precision <- length(intersect(selected_ids, true_ids)) /
      length(unique(selected_ids))
    safe_metric_divide(early_precision, k / n_possible, default = NA_real_)
  }

  act <- signed_epr_one(activation_ids, inhibitory_ids)
  inh <- signed_epr_one(inhibitory_ids, activation_ids)
  metrics <- data.frame(
    Metric = c("EPRActivation", "EPRInhibitory"),
    Value = format_metric_values(c(act, inh))
  )
  list(metrics = metrics)
}

#' @title Calculate Motif Ratios
#' @description Calculates mutual interaction, feedforward-loop and feedback-loop ratios against the ground truth network.
#' These ratios compare local network motif counts in the predicted top-ranked
#' network against the corresponding motif counts in the reference network.
#' @param network_table A data frame of predicted network structure
#' @param ground_truth A data frame of ground truth network
#' @param top_k Number of top-ranked edges to retain. Defaults to the number of ground-truth edges.
#' @param strategy Edge-selection strategy. Use `"beeline"` to match the
#' BEELINE motif evaluator, or `"default"` to use the package's general
#' top-ranked edge selector. Default is `"beeline"`.
#' @return A list containing motif ratios
#' @export
calculate_motif_ratios <- function(
  network_table,
  ground_truth,
  top_k = NULL,
  strategy = "beeline"
) {
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("calculate_motif_ratios() requires the 'igraph' package")
  }
  strategy <- match.arg(strategy, c("beeline", "default"))

  truth <- normalize_ground_truth_edges(ground_truth)
  if (is.null(top_k)) {
    top_k <- nrow(truth)
  }
  pred <- if (identical(strategy, "beeline")) {
    pred_raw <- as.data.frame(network_table, stringsAsFactors = FALSE)
    keep <- intersect(c("regulator", "target", "weight"), names(pred_raw))
    if (length(keep) < 3) {
      stop("network_table must contain columns: regulator, target, weight")
    }
    pred_raw <- pred_raw[, keep, drop = FALSE]
    names(pred_raw) <- c("regulator", "target", "weight")
    pred_raw$weight <- suppressWarnings(as.numeric(pred_raw$weight))
    pred_raw <- pred_raw[
      is.finite(pred_raw$weight) &
        !is.na(pred_raw$regulator) &
        !is.na(pred_raw$target) &
        nzchar(pred_raw$regulator) &
        nzchar(pred_raw$target) &
        pred_raw$regulator != pred_raw$target, ,
      drop = FALSE
    ]
    pred_raw$abs_weight <- abs(pred_raw$weight)
    pred_raw <- pred_raw[
      order(pred_raw$abs_weight, decreasing = TRUE), ,
      drop = FALSE
    ]
    utils::head(pred_raw, top_k)
  } else {
    select_top_ranked_edges(network_table, top_k = top_k)
  }

  motif_counts <- function(edge_df) {
    if (!nrow(edge_df)) {
      return(c(FBL = 0, FFL = 0, MI = 0))
    }
    graph <- igraph::graph_from_data_frame(
      edge_df[, c("regulator", "target"), drop = FALSE],
      directed = TRUE
    )
    graph <- igraph::simplify(
      graph,
      remove.loops = TRUE,
      remove.multiple = TRUE
    )
    edges_mat <- if (igraph::ecount(graph) > 0) {
      igraph::as_edgelist(graph, names = TRUE)
    } else {
      matrix(character(0), ncol = 2)
    }
    edge_ids <- if (length(edges_mat)) {
      apply(edges_mat, 1, paste, collapse = "\r")
    } else {
      character(0)
    }
    reciprocal_pairs <- if (length(edge_ids)) {
      sum(paste(edges_mat[, 2], edges_mat[, 1], sep = "\r") %in% edge_ids) / 2
    } else {
      0
    }

    successors <- split(edges_mat[, 2], edges_mat[, 1])
    for (node in igraph::V(graph)$name) {
      if (is.null(successors[[node]])) {
        successors[[node]] <- character(0)
      } else {
        successors[[node]] <- unique(successors[[node]])
      }
    }

    feedforward_loops <- 0
    for (a in names(successors)) {
      out_a <- successors[[a]]
      if (!length(out_a)) {
        next
      }
      for (b in out_a) {
        out_b <- successors[[b]]
        if (!length(out_b)) {
          next
        }
        for (c in out_b) {
          if (c %in% out_a && c != a && c != b) {
            feedforward_loops <- feedforward_loops + 1
          }
        }
      }
    }

    feedback_loops <- 0
    if (igraph::ecount(graph) > 0) {
      cycles <- igraph::simple_cycles(graph, min = 3, max = 3)
      feedback_loops <- length(cycles)
    }

    c(
      FBL = feedback_loops,
      FFL = feedforward_loops,
      MI = reciprocal_pairs
    )
  }

  pred_counts <- motif_counts(pred)
  truth_counts <- motif_counts(truth)
  ratios <- mapply(
    safe_metric_divide,
    pred_counts,
    truth_counts,
    MoreArgs = list(default = NA_real_)
  )

  list(
    metrics = data.frame(
      Metric = c("FBLRatio", "FFLRatio", "MIRatio"),
      Value = format_metric_values(ratios)
    )
  )
}

#' @title Calculate Path Statistics
#' @description Characterises false-positive predicted edges by their shortest-path distances in the ground truth network.
#' This helps distinguish arbitrary false positives from edges that are close to
#' the true network topology, such as shortcut edges across 2-5 step paths.
#' @param network_table A data frame of predicted network structure
#' @param ground_truth A data frame of ground truth network
#' @param top_k Number of top-ranked edges to retain. Defaults to the number of ground-truth edges.
#' @return A list containing path statistics
#' @export
calculate_path_stats <- function(network_table, ground_truth, top_k = NULL) {
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("calculate_path_stats() requires the 'igraph' package")
  }

  truth <- normalize_ground_truth_edges(ground_truth)
  if (is.null(top_k)) {
    top_k <- nrow(truth)
  }
  pred <- select_top_ranked_edges(
    network_table,
    top_k = top_k,
    drop_zero = TRUE
  )

  ref_graph <- igraph::graph_from_data_frame(
    truth[, c("regulator", "target"), drop = FALSE],
    directed = TRUE
  )
  ref_graph <- igraph::simplify(
    ref_graph,
    remove.loops = TRUE,
    remove.multiple = TRUE
  )
  pred_graph <- if (nrow(pred)) {
    igraph::graph_from_data_frame(
      pred[, c("regulator", "target"), drop = FALSE],
      directed = TRUE
    )
  } else {
    igraph::make_empty_graph(directed = TRUE)
  }
  pred_graph <- igraph::simplify(
    pred_graph,
    remove.loops = TRUE,
    remove.multiple = TRUE
  )

  ref_edges <- if (igraph::ecount(ref_graph)) {
    apply(
      igraph::as_edgelist(ref_graph, names = TRUE),
      1,
      paste,
      collapse = "\r"
    )
  } else {
    character(0)
  }
  pred_edges_mat <- if (igraph::ecount(pred_graph)) {
    igraph::as_edgelist(pred_graph, names = TRUE)
  } else {
    matrix(character(0), ncol = 2)
  }
  pred_edges <- if (length(pred_edges_mat)) {
    apply(pred_edges_mat, 1, paste, collapse = "\r")
  } else {
    character(0)
  }
  false_positive_ids <- setdiff(pred_edges, ref_edges)

  path_counts <- c(`0` = 0, `2` = 0, `3` = 0, `4` = 0, `5` = 0)
  fp_with_path <- 0
  fp_no_path <- 0

  if (length(false_positive_ids)) {
    for (edge_id in false_positive_ids) {
      parts <- strsplit(edge_id, "\r", fixed = TRUE)[[1]]
      if (!all(parts %in% igraph::V(ref_graph)$name)) {
        fp_no_path <- fp_no_path + 1
        next
      }
      distance <- suppressWarnings(igraph::distances(
        ref_graph,
        v = parts[[1]],
        to = parts[[2]],
        mode = "out"
      )[1, 1])
      if (is.infinite(distance)) {
        fp_no_path <- fp_no_path + 1
      } else {
        fp_with_path <- fp_with_path + 1
        key <- as.character(distance)
        if (key %in% names(path_counts)) {
          path_counts[[key]] <- path_counts[[key]] + 1
        }
      }
    }
  }

  metrics <- data.frame(
    Metric = c(
      names(path_counts),
      "numPred",
      "numTP",
      "numFP_withPath",
      "numFP_noPath"
    ),
    Value = c(
      unname(path_counts),
      length(pred_edges),
      length(intersect(pred_edges, ref_edges)),
      fp_with_path,
      fp_no_path
    )
  )
  list(metrics = metrics)
}

#' @title Calculate Jaccard Stability Across Runs
#' @description Calculates the median and dispersion of pairwise top-k Jaccard overlap across multiple predicted networks.
#' This metric measures edge-set stability across repeated runs or perturbations.
#' @param network_tables A list of predicted network tables
#' @param ground_truth A data frame of ground truth network
#' @return A list containing Jaccard stability metrics
#' @export
calculate_stability_jaccard <- function(network_tables, ground_truth) {
  if (!is.list(network_tables) || !length(network_tables)) {
    stop("network_tables must be a non-empty list of predicted networks")
  }

  truth <- normalize_ground_truth_edges(ground_truth)
  top_k <- nrow(truth)
  edge_sets <- lapply(network_tables, function(x) {
    pred <- select_top_ranked_edges(x, top_k = top_k)
    paste(pred$regulator, pred$target, sep = "\r")
  })

  if (length(edge_sets) < 2) {
    scores <- numeric(0)
  } else {
    pair_index <- utils::combn(seq_along(edge_sets), 2)
    scores <- apply(pair_index, 2, function(idx) {
      set_a <- unique(edge_sets[[idx[[1]]]])
      set_b <- unique(edge_sets[[idx[[2]]]])
      union_n <- length(unique(c(set_a, set_b)))
      if (union_n == 0) {
        NA_real_
      } else {
        length(intersect(set_a, set_b)) / union_n
      }
    })
  }

  list(metrics = pairwise_summary(scores, "MedianJaccard", "MADJaccard"))
}

#' @title Calculate Spearman Stability Across Runs
#' @description Calculates the median and dispersion of pairwise Spearman correlations across multiple predicted networks.
#' This metric measures ranking stability across repeated runs or perturbations.
#' @param network_tables A list of predicted network tables
#' @param ground_truth A data frame of ground truth network
#' @return A list containing Spearman stability metrics
#' @export
calculate_stability_spearman <- function(network_tables, ground_truth) {
  if (!is.list(network_tables) || !length(network_tables)) {
    stop("network_tables must be a non-empty list of predicted networks")
  }

  truth <- normalize_ground_truth_edges(ground_truth)
  gene_set <- sort(unique(c(truth$regulator, truth$target)))
  if (length(gene_set) < 2) {
    correlations <- numeric(0)
  } else {
    universe <- expand.grid(
      regulator = gene_set,
      target = gene_set,
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )
    universe <- universe[universe$regulator != universe$target, , drop = FALSE]
    universe_ids <- paste(universe$regulator, universe$target, sep = "\r")

    weight_matrix <- lapply(network_tables, function(x) {
      pred <- normalize_signed_predicted_edges(x)
      pred_ids <- paste(pred$regulator, pred$target, sep = "\r")
      vec <- pred$weight[match(universe_ids, pred_ids)]
      vec[is.na(vec)] <- 0
      vec
    })
    weight_matrix <- Filter(function(x) !all(x == x[[1]]), weight_matrix)

    if (length(weight_matrix) < 2) {
      correlations <- numeric(0)
    } else {
      matrix_data <- do.call(cbind, weight_matrix)
      corr_matrix <- stats::cor(matrix_data, method = "spearman")
      correlations <- corr_matrix[upper.tri(corr_matrix)]
      correlations <- correlations[is.finite(correlations)]
    }
  }

  list(
    metrics = pairwise_summary(correlations, "MedianSpearman", "MADSpearman")
  )
}

plot_all_metrics <- function(plot_metrics, line_color = "#1563cc") {
  plot_metrics$Metric <- factor(
    plot_metrics$Metric,
    levels = c(
      "AUROC",
      "AUPRC",
      "EPR",
      "Precision",
      "Recall",
      "F1",
      "Accuracy",
      "JI"
    )
  )
  ggplot2::ggplot(
    plot_metrics,
    ggplot2::aes(x = Metric, y = Value, fill = Metric)
  ) +
    ggplot2::geom_col(color = line_color, linewidth = 0.3) +
    ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%.3f", Value)),
      vjust = -0.2,
      size = 3
    ) +
    ggplot2::scale_fill_manual(
      values = rep(line_color, length(unique(plot_metrics$Metric)))
    ) +
    ggplot2::coord_cartesian(
      ylim = c(0, max(plot_metrics$Value, na.rm = TRUE) * 1.1)
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "none")
}

#' @title Plot Edges Comparison
#' @description Creates a scatter plot comparing predicted, ground truth, and
#'   overlapping edges between two gene regulatory networks.
#' @param network_table A data frame of predicted network structure with
#'   `regulator` and `target` columns.
#' @param ground_truth A data frame of ground truth network with `regulator`
#'   and `target` columns.
#' @param color_pattern A named list of colors for predicted, ground_truth,
#'   overlap, and total edges.
#' @return A ggplot object visualizing edge overlap between the two networks.
#' @export
plot_edges_comparison <- function(
  network_table,
  ground_truth,
  color_pattern = list(
    predicted = "gray",
    ground_truth = "#bb141a",
    overlap = "#1966ad",
    total = "#6C757D"
  )
) {
  all_genes <- sort(
    unique(
      c(
        network_table$regulator,
        network_table$target,
        ground_truth$regulator,
        ground_truth$target
      )
    )
  )

  pred_edges <- network_table[, c("regulator", "target")]
  true_edges <- ground_truth[, c("regulator", "target")]

  pred_edge_ids <- paste(
    pred_edges$regulator,
    pred_edges$target,
    sep = "-"
  )
  true_edge_ids <- paste(
    true_edges$regulator,
    true_edges$target,
    sep = "-"
  )
  overlap_edges <- intersect(
    pred_edge_ids,
    true_edge_ids
  )

  categories <- c(
    "Predicted Only" = "Predicted",
    "Ground Truth Only" = "Ground Truth",
    "Overlapping" = "Overlap"
  )

  edge_plot_data <- rbind(
    data.frame(
      regulator = pred_edges$regulator,
      target = pred_edges$target,
      type = "Predicted Only",
      stringsAsFactors = FALSE
    ),
    data.frame(
      regulator = true_edges$regulator,
      target = true_edges$target,
      type = "Ground Truth Only",
      stringsAsFactors = FALSE
    )
  )

  if (length(overlap_edges) > 0) {
    overlap_data <- data.frame(
      regulator = sub("-.*", "", overlap_edges),
      target = sub(".*-", "", overlap_edges),
      type = "Overlapping",
      stringsAsFactors = FALSE
    )
    edge_plot_data <- rbind(edge_plot_data, overlap_data)
  }

  edge_plot_data$type <- factor(
    edge_plot_data$type,
    levels = c(
      "Predicted Only",
      "Ground Truth Only",
      "Overlapping"
    )
  )

  edge_plot_data$regulator <- factor(
    edge_plot_data$regulator,
    levels = all_genes
  )
  edge_plot_data$target <- factor(
    edge_plot_data$target,
    levels = all_genes
  )

  network_plot <- ggplot() +
    geom_tile(
      data = edge_plot_data,
      aes(y = regulator, x = target, fill = type),
      color = "white",
      width = 0.9,
      height = 0.9
    ) +
    scale_fill_manual(
      values = c(
        "Predicted Only" = color_pattern$predicted,
        "Ground Truth Only" = color_pattern$ground_truth,
        "Overlapping" = color_pattern$overlap
      ),
      labels = categories
    ) +
    theme_bw() +
    labs(
      title = "Network Edge Comparison",
      y = "Regulator",
      x = "Target"
    ) +
    theme(
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        vjust = 1
      ),
      legend.position = "bottom",
      legend.title = element_blank(),
      panel.grid = element_line(color = "gray95"),
      plot.margin = margin(b = 20)
    ) +
    coord_fixed()

  intersection_size <- length(overlap_edges)
  predicted_only <- length(setdiff(pred_edge_ids, true_edge_ids))
  ground_truth_only <- length(setdiff(true_edge_ids, pred_edge_ids))
  total_predicted <- length(pred_edge_ids)
  total_ground_truth <- length(true_edge_ids)
  ordered_categories <- c(
    "Total Predicted",
    "Total Ground Truth",
    "Predicted Only",
    "Ground Truth Only",
    "Overlapping"
  )
  edge_details <- data.frame(
    Category = factor(
      ordered_categories,
      levels = ordered_categories
    ),
    count = c(
      total_predicted,
      total_ground_truth,
      predicted_only,
      ground_truth_only,
      intersection_size
    ),
    type = factor(
      c(
        "Total",
        "Total",
        "Predicted",
        "Ground Truth",
        "Overlap"
      ),
      levels = c(
        "Total",
        "Predicted",
        "Ground Truth",
        "Overlap"
      )
    )
  )

  edge_stats_plot <- ggplot(
    edge_details,
    aes(x = Category, y = count, fill = type)
  ) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_text(aes(label = count), vjust = -0.3) +
    scale_fill_manual(
      values = c(
        "Total" = color_pattern$total,
        "Predicted" = color_pattern$predicted,
        "Ground Truth" = color_pattern$ground_truth,
        "Overlap" = color_pattern$overlap
      )
    ) +
    theme_bw() +
    labs(
      title = "Edge Distribution",
      x = "",
      y = "Count"
    ) +
    theme(
      axis.text.x = element_text(angle = 30, hjust = 1),
      legend.position = "none",
      legend.title = element_blank()
    ) +
    scale_y_continuous(
      expand = expansion(mult = c(0, 0.1)),
      limits = function(x) c(0, max(x) * 1.05)
    )

  final_plot <- patchwork::wrap_plots(
    list(network_plot, edge_stats_plot),
    ncol = 2
  ) +
    patchwork::plot_annotation(tag_levels = "A")

  final_plot
}
