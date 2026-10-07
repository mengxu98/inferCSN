#include <Rcpp.h>
#include <algorithm>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>
using namespace Rcpp;

bool geneCompare(const std::string &a, const std::string &b)
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
}

//' @title Filter and sort a network matrix
//'
//' @param network_matrix Network weight matrix.
//' @param regulators,targets Nodes to include.
//' @return A filtered and sorted matrix.
//' @export
// [[Rcpp::export]]
NumericMatrix
filter_sort_matrix(NumericMatrix network_matrix,
                   Nullable<CharacterVector> regulators = R_NilValue,
                   Nullable<CharacterVector> targets = R_NilValue)
{
  for (R_xlen_t i = 0; i < network_matrix.length(); i++)
  {
    if (R_IsNA(network_matrix[i]))
    {
      network_matrix[i] = 0;
    }
  }

  CharacterVector curr_regulators = rownames(network_matrix);
  CharacterVector curr_targets = colnames(network_matrix);

  const auto select_indices = [](CharacterVector current,
                                 Nullable<CharacterVector> requested) {
    std::vector<std::string> selected;
    std::vector<std::string> allowed;
    if (requested.isNotNull()) allowed = as<std::vector<std::string>>(CharacterVector(requested));
    std::unordered_map<std::string, int> positions;
    for (int i = 0; i < current.size(); ++i) {
      const std::string name = as<std::string>(current[i]);
      positions[name] = i;
      if (requested.isNull() || std::find(allowed.begin(), allowed.end(), name) != allowed.end()) {
        selected.push_back(name);
      }
    }
    std::sort(selected.begin(), selected.end(), geneCompare);
    IntegerVector indices(selected.size());
    for (size_t i = 0; i < selected.size(); ++i) indices[i] = positions[selected[i]];
    CharacterVector names = wrap(selected);
    return std::make_pair(indices, names);
  };
  const auto selected_rows = select_indices(curr_regulators, regulators);
  const auto selected_columns = select_indices(curr_targets, targets);
  IntegerVector rows = selected_rows.first;
  IntegerVector columns = selected_columns.first;
  NumericMatrix result(rows.size(), columns.size());
  for (int i = 0; i < rows.size(); ++i)
  {
    for (int j = 0; j < columns.size(); ++j)
    {
      result(i, j) = network_matrix(rows[i], columns[j]);
    }
  }
  rownames(result) = selected_rows.second;
  colnames(result) = selected_columns.second;

  return result;
}
