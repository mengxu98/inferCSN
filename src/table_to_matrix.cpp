#include <Rcpp.h>
#include <algorithm>
#include <unordered_map>

using namespace Rcpp;

//' @title Switch network table to matrix
//'
//' @param network_table The weight data table of network.
//' @param regulators Regulators list.
//' @param targets Targets list.
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
//' table_to_matrix(
//'   network_table,
//'   regulators = c("g1", "g2"),
//'   targets = c("g3", "g4")
//' )
// [[Rcpp::export]]
NumericMatrix table_to_matrix(DataFrame network_table,
                              Nullable<CharacterVector> regulators = R_NilValue,
                              Nullable<CharacterVector> targets = R_NilValue)
{
  CharacterVector table_regulators = network_table["regulator"];
  CharacterVector table_targets = network_table["target"];
  NumericVector weight = network_table["weight"];

  // Handle optional regulators and targets parameters
  CharacterVector filter_regulators;
  CharacterVector filter_targets;

  if (regulators.isNotNull())
  {
    CharacterVector reg(regulators);
    filter_regulators = intersect(unique(table_regulators), reg);
  }
  else
  {
    filter_regulators = unique(table_regulators);
  }

  if (targets.isNotNull())
  {
    CharacterVector tar(targets);
    filter_targets = intersect(unique(table_targets), tar);
  }
  else
  {
    filter_targets = unique(table_targets);
  }

  // Convert to standard strings for sorting
  std::vector<std::string> reg_strings;
  std::vector<std::string> tar_strings;

  for (size_t i = 0; i < filter_regulators.length(); i++)
  {
    reg_strings.push_back(Rcpp::as<std::string>(filter_regulators[i]));
  }
  for (size_t i = 0; i < filter_targets.length(); i++)
  {
    tar_strings.push_back(Rcpp::as<std::string>(filter_targets[i]));
  }

  // Custom comparison function for gene names (e.g., g1, g2, g10)
  auto geneCompare = [](const std::string &a, const std::string &b)
  {
    size_t na = a.find_first_of("0123456789");
    size_t nb = b.find_first_of("0123456789");

    if (na != std::string::npos && nb != std::string::npos)
    {
      std::string prefix_a = a.substr(0, na);
      std::string prefix_b = b.substr(0, nb);
      if (prefix_a == prefix_b)
      {
        return std::stoi(a.substr(na)) < std::stoi(b.substr(nb));
      }
    }
    return a < b;
  };

  // Sort strings
  std::sort(reg_strings.begin(), reg_strings.end(), geneCompare);
  std::sort(tar_strings.begin(), tar_strings.end(), geneCompare);

  // Create maps for lookups
  std::unordered_map<std::string, int> regulator_indices;
  std::unordered_map<std::string, int> target_indices;

  for (size_t i = 0; i < reg_strings.size(); ++i)
  {
    regulator_indices[reg_strings[i]] = i;
  }
  for (size_t i = 0; i < tar_strings.size(); ++i)
  {
    target_indices[tar_strings[i]] = i;
  }

  // Create and initialize matrix with zeros
  NumericMatrix weight_matrix(reg_strings.size(), tar_strings.size());
  std::fill(weight_matrix.begin(), weight_matrix.end(), 0.0);

  // Convert back to CharacterVector for rownames/colnames
  CharacterVector sorted_regulators(reg_strings.size());
  CharacterVector sorted_targets(tar_strings.size());
  for (size_t i = 0; i < reg_strings.size(); i++)
  {
    sorted_regulators[i] = reg_strings[i];
  }
  for (size_t i = 0; i < tar_strings.size(); i++)
  {
    sorted_targets[i] = tar_strings[i];
  }

  rownames(weight_matrix) = sorted_regulators;
  colnames(weight_matrix) = sorted_targets;

  // Fill matrix only with filtered and valid entries
  for (size_t i = 0; i < network_table.nrows(); ++i)
  {
    std::string reg = Rcpp::as<std::string>(table_regulators[i]);
    std::string tar = Rcpp::as<std::string>(table_targets[i]);

    auto reg_it = regulator_indices.find(reg);
    auto tar_it = target_indices.find(tar);

    // Skip if regulator or target is not in filtered set
    if (reg_it != regulator_indices.end() && tar_it != target_indices.end())
    {
      weight_matrix(reg_it->second, tar_it->second) = weight[i];
    }
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
data("example_matrix")
network_table <- inferCSN(
  example_matrix,
  verbose = TRUE
)

weight_matrix_c1 <- table_to_matrix(network_table)
weight_matrix_r1 <- table_to_matrix_v1(network_table)
weight_matrix_r2 <- table_to_matrix_v2(network_table)

identical(weight_matrix, weight_matrix_r1)
identical(weight_matrix, weight_matrix_r2)

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
