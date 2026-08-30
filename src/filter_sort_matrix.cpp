#include <Rcpp.h>
#include <algorithm>
#include <string>
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

  std::vector<std::string> filtered_regulators;
  if (regulators.isNotNull())
  {
    CharacterVector reg(regulators);
    for (R_xlen_t i = 0; i < curr_regulators.length(); i++)
    {
      std::string curr_reg = as<std::string>(curr_regulators[i]);
      for (R_xlen_t j = 0; j < reg.length(); j++)
      {
        if (curr_reg == as<std::string>(reg[j]))
        {
          filtered_regulators.push_back(curr_reg);
          break;
        }
      }
    }
  }
  else
  {
    for (R_xlen_t i = 0; i < curr_regulators.length(); i++)
    {
      filtered_regulators.push_back(as<std::string>(curr_regulators[i]));
    }
  }

  std::vector<std::string> filtered_targets;
  if (targets.isNotNull())
  {
    CharacterVector tar(targets);
    for (R_xlen_t i = 0; i < curr_targets.length(); i++)
    {
      std::string curr_tar = as<std::string>(curr_targets[i]);
      for (R_xlen_t j = 0; j < tar.length(); j++)
      {
        if (curr_tar == as<std::string>(tar[j]))
        {
          filtered_targets.push_back(curr_tar);
          break;
        }
      }
    }
  }
  else
  {
    for (R_xlen_t i = 0; i < curr_targets.length(); i++)
    {
      filtered_targets.push_back(as<std::string>(curr_targets[i]));
    }
  }

  std::sort(filtered_regulators.begin(), filtered_regulators.end(),
            geneCompare);
  std::sort(filtered_targets.begin(), filtered_targets.end(), geneCompare);

  NumericMatrix result(filtered_regulators.size(), filtered_targets.size());

  std::unordered_map<std::string, int> old_reg_indices;
  std::unordered_map<std::string, int> old_tar_indices;

  for (R_xlen_t i = 0; i < curr_regulators.length(); i++)
  {
    old_reg_indices[as<std::string>(curr_regulators[i])] = i;
  }
  for (R_xlen_t i = 0; i < curr_targets.length(); i++)
  {
    old_tar_indices[as<std::string>(curr_targets[i])] = i;
  }

  for (size_t i = 0; i < filtered_regulators.size(); i++)
  {
    for (size_t j = 0; j < filtered_targets.size(); j++)
    {
      int old_row = old_reg_indices[filtered_regulators[i]];
      int old_col = old_tar_indices[filtered_targets[j]];
      result(i, j) = network_matrix(old_row, old_col);
    }
  }

  CharacterVector new_regulators(filtered_regulators.size());
  CharacterVector new_targets(filtered_targets.size());

  for (size_t i = 0; i < filtered_regulators.size(); i++)
  {
    new_regulators[i] = filtered_regulators[i];
  }
  for (size_t i = 0; i < filtered_targets.size(); i++)
  {
    new_targets[i] = filtered_targets[i];
  }

  rownames(result) = new_regulators;
  colnames(result) = new_targets;

  return result;
}
