#include <Rcpp.h>
#include <algorithm>
#include <unordered_map>
#include <vector>

using namespace Rcpp;

//' @title Switch network table to matrix
//'
//' @param network_table The weight data table of network.
//' @param regulators Regulators list.
//' @param targets Targets list.
//' @param threshold The threshold for filtering weights based on absolute values, defaults to 0.
//'
//' @return Weight matrix
//' @export
//'
//' @examples
//' data("example_matrix")
//' network_table <- inferCSN(example_matrix)
//' head(network_table)
//'
//' table_to_matrix(network_table)[1:6, 1:6]
//'
//' table_to_matrix(network_table, threshold = 0.8)
//'
//' table_to_matrix(
//'   network_table,
//'   regulators = c("g1", "g2"),
//'   targets = c("g3", "g4")
//' )
// [[Rcpp::export]]
NumericMatrix table_to_matrix(DataFrame network_table,
                              Nullable<CharacterVector> regulators = R_NilValue,
                              Nullable<CharacterVector> targets = R_NilValue,
                              double threshold = 0.0) {
  CharacterVector table_regulators = network_table["regulator"];
  CharacterVector table_targets = network_table["target"];
  NumericVector weight = network_table["weight"];

  std::unordered_map<std::string, size_t> unique_regs;
  std::unordered_map<std::string, size_t> unique_tars;
  std::vector<std::pair<size_t, size_t>> valid_indices;
  std::vector<double> valid_weights;
  valid_indices.reserve(weight.length());
  valid_weights.reserve(weight.length());

  for (int i = 0; i < weight.length(); i++) {
    if (std::abs(weight[i]) >= threshold) {

      String reg_str = table_regulators[i];
      String tar_str = table_targets[i];
      std::string reg = reg_str;
      std::string tar = tar_str;

      if (unique_regs.find(reg) == unique_regs.end()) {
        unique_regs[reg] = unique_regs.size();
      }
      if (unique_tars.find(tar) == unique_tars.end()) {
        unique_tars[tar] = unique_tars.size();
      }

      valid_indices.emplace_back(unique_regs[reg], unique_tars[tar]);
      valid_weights.push_back(weight[i]);
    }
  }

  std::vector<std::string> reg_strings;
  std::vector<std::string> tar_strings;
  reg_strings.reserve(unique_regs.size());
  tar_strings.reserve(unique_tars.size());

  for (const auto &pair : unique_regs) {
    reg_strings.push_back(pair.first);
  }
  for (const auto &pair : unique_tars) {
    tar_strings.push_back(pair.first);
  }

  auto geneCompare = [](const std::string &a, const std::string &b) {
    size_t na = a.find_first_of("0123456789");
    size_t nb = b.find_first_of("0123456789");

    if (na != std::string::npos && nb != std::string::npos) {
      std::string prefix_a = a.substr(0, na);
      std::string prefix_b = b.substr(0, nb);
      if (prefix_a == prefix_b) {
        return std::stoi(a.substr(na)) < std::stoi(b.substr(nb));
      }
    }
    return a < b;
  };

  std::sort(reg_strings.begin(), reg_strings.end(), geneCompare);
  std::sort(tar_strings.begin(), tar_strings.end(), geneCompare);

  if (regulators.isNotNull()) {
    CharacterVector reg_filter(regulators);
    std::vector<std::string> filtered_regs;
    filtered_regs.reserve(reg_strings.size());

    for (int i = 0; i < reg_filter.length(); i++) {
      String reg_str = reg_filter[i];
      std::string reg = reg_str;
      if (unique_regs.find(reg) != unique_regs.end()) {
        filtered_regs.push_back(reg);
      }
    }
    reg_strings = std::move(filtered_regs);
  }

  if (targets.isNotNull()) {
    CharacterVector tar_filter(targets);
    std::vector<std::string> filtered_tars;
    filtered_tars.reserve(tar_strings.size());

    for (int i = 0; i < tar_filter.length(); i++) {
      String tar_str = tar_filter[i];
      std::string tar = tar_str;
      if (unique_tars.find(tar) != unique_tars.end()) {
        filtered_tars.push_back(tar);
      }
    }
    tar_strings = std::move(filtered_tars);
  }

  NumericMatrix weight_matrix(reg_strings.size(), tar_strings.size());
  std::fill(weight_matrix.begin(), weight_matrix.end(), 0.0);

  CharacterVector row_names(reg_strings.size());
  CharacterVector col_names(tar_strings.size());

  std::unordered_map<std::string, size_t> reg_indices;
  std::unordered_map<std::string, size_t> tar_indices;

  for (size_t i = 0; i < reg_strings.size(); ++i) {
    row_names[i] = reg_strings[i];
    reg_indices[reg_strings[i]] = i;
  }
  for (size_t i = 0; i < tar_strings.size(); ++i) {
    col_names[i] = tar_strings[i];
    tar_indices[tar_strings[i]] = i;
  }

  rownames(weight_matrix) = row_names;
  colnames(weight_matrix) = col_names;

  for (size_t i = 0; i < valid_indices.size(); ++i) {
    const auto &[reg_idx, tar_idx] = valid_indices[i];
    weight_matrix(reg_idx, tar_idx) = valid_weights[i];
  }

  return weight_matrix;
}

/*
table_to_matrix_v1 <- function(
    network_table,
    regulators = NULL,
    targets = NULL) {
  network_table <- network_format(
    network_table,
    regulators = regulators,
    targets = targets,
    abs_weight = FALSE
  )
  unique_regulators <- gtools::mixedsort(
    unique(network_table$regulator)
  )
  unique_targets <- gtools::mixedsort(
    unique(network_table$target)
  )

  weight_matrix <- matrix(
    0,
    nrow = length(unique_regulators),
    ncol = length(unique_targets),
    dimnames = list(unique_regulators, unique_targets)
  )

  for (i in 1:nrow(network_table)) {
    weight_matrix[network_table$regulator[i], network_table$target[i]] <-
      network_table$weight[i]
  }

  return(weight_matrix)
}

table_to_matrix_v2 <- function(
    network_table,
    regulators = NULL,
    targets = NULL) {
  network_table <- network_format(
    network_table,
    regulators = regulators,
    targets = targets,
    abs_weight = FALSE
  )
  unique_regulators <- gtools::mixedsort(
    unique(network_table$regulator)
  )
  unique_targets <- gtools::mixedsort(
    unique(network_table$target)
  )
  weight_matrix <- suppressMessages(
    reshape2::acast(
      network_table,
      regulator ~ target,
      fill = 0
    )[unique_regulators, unique_targets]
  )

  return(weight_matrix)
}

devtools::document()
network_table <- data.frame(
  "regulator" = paste0("g", 1:10000),
  "target" = paste0("g", 5001:15000),
  "weight" = runif(10000)
)

weight_matrix_c1 <- table_to_matrix(network_table)
weight_matrix_r1 <- table_to_matrix_v1(network_table)
weight_matrix_r2 <- table_to_matrix_v2(network_table)

identical(weight_matrix_c1, weight_matrix_r1)
identical(weight_matrix_c1, weight_matrix_r2)

weight_matrix_c2 <- table_to_matrix(
  network_table,
  regulators = c("g1", "g2", "g3"),
  targets = c("g4", "g5", "g6")
)
weight_matrix_r3 <- table_to_matrix_v1(
  network_table,
  regulators = c("g1", "g2", "g3"),
  targets = c("g4", "g5", "g6")
)
weight_matrix_r4 <- table_to_matrix_v2(
  network_table,
  regulators = c("g1", "g2", "g3"),
  targets = c("g4", "g5", "g6")
)
identical(weight_matrix_c2, weight_matrix_r3)
identical(weight_matrix_c2, weight_matrix_r4)

if (!requireNamespace("bench", quietly = TRUE)) {
  install.packages("bench")
}
if (!requireNamespace("reshape2", quietly = TRUE)) {
  install.packages("reshape2")
}

bench::mark(
  table_to_matrix(network_table),
  table_to_matrix_v1(network_table),
  table_to_matrix_v2(network_table)
)
# # A tibble: 3 × 13
#   expression               min   median `itr/sec` mem_alloc `gc/sec` n_itr  n_gc
#   <bch:expr>          <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl> <int> <dbl>
# 1 table_to_matrix(ne… 222.04ms 349.54ms     2.86    763.1MB    1.43      2     1
# 2 table_to_matrix_v1…     1.3s     1.3s     0.767    2.25GB    0.767     1     1
# 3 table_to_matrix_v2…    2.09s    2.09s     0.479    4.12GB    0.959     1     2
 */
