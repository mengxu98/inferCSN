#include "network_format.h"
#include <Rcpp.h>
#include <algorithm>

using namespace Rcpp;

//' @title Switch matrix to network table
//'
//' @param network_matrix The matrix of network weight.
//' @param regulators Character vector of regulator names to filter by.
//' @param targets Character vector of target names to filter by.
//' @param threshold The threshold for filtering weights based on absolute values, defaults to 0.
//'
//' @return Network table
//' @export
//'
//' @examples
//' data("example_matrix")
//' network_table <- inferCSN(example_matrix)
//' network_matrix <- table_to_matrix(network_table)
//' network_table_new <- matrix_to_table(network_matrix)
//' head(network_table)
//' head(network_table_new)
//' identical(
//'   network_table,
//'   network_table_new
//' )
//'
//' matrix_to_table(
//'   network_matrix,
//'   threshold = 0.8
//' )
//'
//' matrix_to_table(
//'   network_matrix,
//'   regulators = c("g1", "g2"),
//'   targets = c("g3", "g4")
//' )
//'
// [[Rcpp::export]]
DataFrame matrix_to_table(NumericMatrix network_matrix,
                          Nullable<CharacterVector> regulators = R_NilValue,
                          Nullable<CharacterVector> targets = R_NilValue,
                          double threshold = 0.0) {
  int nrow = network_matrix.nrow();
  int ncol = network_matrix.ncol();

  CharacterVector row_names = rownames(network_matrix);
  CharacterVector col_names = colnames(network_matrix);

  if (row_names.length() == 0 || col_names.length() == 0) {
    stop("Input matrix must have both row and column names");
  }

  CharacterVector reg_filter;
  CharacterVector tar_filter;
  bool use_reg_filter = false;
  bool use_tar_filter = false;

  if (regulators.isNotNull()) {
    reg_filter = as<CharacterVector>(regulators);
    use_reg_filter = true;
  }
  if (targets.isNotNull()) {
    tar_filter = as<CharacterVector>(targets);
    use_tar_filter = true;
  }

  // First count non-zero elements that pass threshold using absolute values
  int valid_count = 0;
  for (int i = 0; i < nrow; ++i) {
    if (use_reg_filter && !std::any_of(reg_filter.begin(), reg_filter.end(),
                                       [&row_names, i](const String &reg) {
                                         return reg == row_names[i];
                                       })) {
      continue;
    }

    for (int j = 0; j < ncol; ++j) {
      if (use_tar_filter && !std::any_of(tar_filter.begin(), tar_filter.end(),
                                         [&col_names, j](const String &tar) {
                                           return tar == col_names[j];
                                         })) {
        continue;
      }

      double weight = network_matrix(i, j);
      if (weight != 0 && std::abs(weight) >= threshold) {
        valid_count++;
      }
    }
  }

  // Pre-allocate vectors with exact size needed
  CharacterVector regulators_out(valid_count);
  CharacterVector targets_out(valid_count);
  NumericVector weights(valid_count);
  IntegerVector indices(valid_count);

  // Fill vectors using absolute values for threshold comparison
  int idx = 0;
  for (int i = 0; i < nrow; ++i) {
    if (use_reg_filter && !std::any_of(reg_filter.begin(), reg_filter.end(),
                                       [&row_names, i](const String &reg) {
                                         return reg == row_names[i];
                                       })) {
      continue;
    }

    for (int j = 0; j < ncol; ++j) {
      if (use_tar_filter && !std::any_of(tar_filter.begin(), tar_filter.end(),
                                         [&col_names, j](const String &tar) {
                                           return tar == col_names[j];
                                         })) {
        continue;
      }

      double weight = network_matrix(i, j);
      if (weight != 0 && std::abs(weight) >= threshold) {
        regulators_out[idx] = row_names[i];
        targets_out[idx] = col_names[j];
        weights[idx] = weight;
        indices[idx] = idx;
        idx++;
      }
    }
  }

  std::sort(indices.begin(), indices.end(), [&weights](int i1, int i2) {
    return std::abs(weights[i1]) > std::abs(weights[i2]);
  });

  CharacterVector sorted_regulators(valid_count);
  CharacterVector sorted_targets(valid_count);
  NumericVector sorted_weights(valid_count);

  for (int i = 0; i < valid_count; i++) {
    sorted_regulators[i] = regulators_out[indices[i]];
    sorted_targets[i] = targets_out[indices[i]];
    sorted_weights[i] = weights[indices[i]];
  }

  DataFrame intermediate =
      DataFrame::create(Rcpp::Named("regulator") = sorted_regulators,
                        Rcpp::Named("target") = sorted_targets,
                        Rcpp::Named("weight") = sorted_weights);

  return network_format(intermediate, regulators, targets, false);
}

/*
matrix_to_table_R <- function(
    weightMatrix,
    reportMax = NULL,
    threshold = 0) {
    if (!is.numeric(threshold)) {
        stop("threshold must be a number.")
    }
    regulatorsInTargets <- rownames(weightMatrix)[rownames(weightMatrix) %in%
colnames(weightMatrix)] if (length(regulatorsInTargets) == 1) {
        weightMatrix[regulatorsInTargets, regulatorsInTargets] <- NA
    }
    if (length(regulatorsInTargets) > 1) {
        diag(weightMatrix[regulatorsInTargets, regulatorsInTargets]) <- NA
    }
    linkList <- reshape2::melt(weightMatrix, na.rm = TRUE)
    colnames(linkList) <- c("regulator", "target", "weight")
    linkList <- linkList[linkList$weight >= threshold, ]
    linkList <- linkList[order(linkList$weight, decreasing = TRUE), ]
    if (!is.null(reportMax)) {
        linkList <- linkList[1:min(nrow(linkList), reportMax), ]
    }
    rownames(linkList) <- NULL
    uniquePairs <- nrow(unique(linkList[, c(
        "regulator",
        "target"
    )]))
    if (uniquePairs < nrow(linkList)) {
        warning("There might be duplicated regulator-target (gene id/name)
pairs.")
    }
    return(linkList)
}

devtools::document()
network_table <- data.frame(
    "regulator" = paste0("g", 1:5000),
    "target" = paste0("g", 2501:7500),
    "weight" = runif(5000)
)

network_matrix <- table_to_matrix(network_table)

network_matrix_C <- matrix_to_table(network_matrix) |>
    network_format(abs_weight = FALSE)

network_matrix_R <- matrix_to_table_R(network_matrix) |>
    network_format(abs_weight = FALSE)

identical(network_matrix_C, network_matrix_R)

bench::mark(
    network_matrix_C <- matrix_to_table(network_matrix) |>
        network_format(abs_weight = FALSE),
    network_matrix_R <- matrix_to_table_R(network_matrix) |>
        network_format(abs_weight = FALSE)
)
# # A tibble: 2 × 13
#   expression     min  median `itr/sec` mem_alloc `gc/sec` n_itr  n_gc   total_time 
#   <bch:expr> <bch:t> <bch:t>     <dbl> <bch:byt>    <dbl> <int> <dbl>   <bch:tm> 
# 1 network_m… 568.2ms 568.2ms    1.76      1.05MB    0         1       0 568.2ms 
# 2 network_m…   33.9s   33.9s    0.0295    9.08GB    0.354     1       1 233.9s
 */
