#include <Rcpp.h>
#include <algorithm>
#include <string>
using namespace Rcpp;

// Custom comparison function for gene names (e.g., g1, g2, g10)
bool geneCompare(const std::string &a, const std::string &b) {
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
}

//' @title Filter and sort matrix
//'
//' @param network_matrix The matrix of network weight.
//' @param regulators Regulators list.
//' @param targets Targets list.
//'
//' @return Filtered and sorted matrix
//' @export
//'
//' @examples
//' data("example_matrix")
//' network_table <- inferCSN(example_matrix)
//' network_matrix <- table_to_matrix(network_table)
//' filter_sort_matrix(network_matrix)[1:6, 1:6]
//'
//' filter_sort_matrix(
//'   network_matrix,
//'   regulators = c("g1", "g2"),
//'   targets = c("g3", "g4")
//' )
// [[Rcpp::export]]
NumericMatrix
filter_sort_matrix(NumericMatrix network_matrix,
                   Nullable<CharacterVector> regulators = R_NilValue,
                   Nullable<CharacterVector> targets = R_NilValue) {
  // Replace NA with 0
  for (size_t i = 0; i < network_matrix.length(); i++) {
    if (R_IsNA(network_matrix[i])) {
      network_matrix[i] = 0;
    }
  }

  // Get current row and column names
  CharacterVector curr_regulators = rownames(network_matrix);
  CharacterVector curr_targets = colnames(network_matrix);

  // Handle regulators filtering
  std::vector<std::string> filtered_regulators;
  if (regulators.isNotNull()) {
    CharacterVector reg(regulators);
    // Get intersection
    for (size_t i = 0; i < curr_regulators.length(); i++) {
      std::string curr_reg = as<std::string>(curr_regulators[i]);
      for (size_t j = 0; j < reg.length(); j++) {
        if (curr_reg == as<std::string>(reg[j])) {
          filtered_regulators.push_back(curr_reg);
          break;
        }
      }
    }
  } else {
    for (size_t i = 0; i < curr_regulators.length(); i++) {
      filtered_regulators.push_back(as<std::string>(curr_regulators[i]));
    }
  }

  // Handle targets filtering
  std::vector<std::string> filtered_targets;
  if (targets.isNotNull()) {
    CharacterVector tar(targets);
    // Get intersection
    for (size_t i = 0; i < curr_targets.length(); i++) {
      std::string curr_tar = as<std::string>(curr_targets[i]);
      for (size_t j = 0; j < tar.length(); j++) {
        if (curr_tar == as<std::string>(tar[j])) {
          filtered_targets.push_back(curr_tar);
          break;
        }
      }
    }
  } else {
    for (size_t i = 0; i < curr_targets.length(); i++) {
      filtered_targets.push_back(as<std::string>(curr_targets[i]));
    }
  }

  // Sort gene names
  std::sort(filtered_regulators.begin(), filtered_regulators.end(),
            geneCompare);
  std::sort(filtered_targets.begin(), filtered_targets.end(), geneCompare);

  // Create new matrix with filtered dimensions
  NumericMatrix result(filtered_regulators.size(), filtered_targets.size());

  // Create maps for quick lookups
  std::unordered_map<std::string, int> old_reg_indices;
  std::unordered_map<std::string, int> old_tar_indices;

  for (size_t i = 0; i < curr_regulators.length(); i++) {
    old_reg_indices[as<std::string>(curr_regulators[i])] = i;
  }
  for (size_t i = 0; i < curr_targets.length(); i++) {
    old_tar_indices[as<std::string>(curr_targets[i])] = i;
  }

  // Fill new matrix
  for (size_t i = 0; i < filtered_regulators.size(); i++) {
    for (size_t j = 0; j < filtered_targets.size(); j++) {
      int old_row = old_reg_indices[filtered_regulators[i]];
      int old_col = old_tar_indices[filtered_targets[j]];
      result(i, j) = network_matrix(old_row, old_col);
    }
  }

  // Set row and column names
  CharacterVector new_regulators(filtered_regulators.size());
  CharacterVector new_targets(filtered_targets.size());

  for (size_t i = 0; i < filtered_regulators.size(); i++) {
    new_regulators[i] = filtered_regulators[i];
  }
  for (size_t i = 0; i < filtered_targets.size(); i++) {
    new_targets[i] = filtered_targets[i];
  }

  rownames(result) = new_regulators;
  colnames(result) = new_targets;

  return result;
}

/*
filter_sort_matrix <- function(
    network_matrix,
    regulators = NULL,
    targets = NULL) {
  network_matrix[is.na(network_matrix)] <- 0
  if (!is.null(regulators)) {
    regulators <- intersect(rownames(network_matrix), regulators)
  } else {
    regulators <- rownames(network_matrix)
  }
  if (!is.null(targets)) {
    targets <- intersect(colnames(network_matrix), targets)
  } else {
    targets <- colnames(network_matrix)
  }

  unique_regulators <- gtools::mixedsort(unique(regulators))
  unique_targets <- gtools::mixedsort(unique(targets))
  network_matrix <- network_matrix[unique_regulators, unique_targets]

  return(network_matrix)
}
*/
