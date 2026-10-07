#include <Rcpp.h>
#include <algorithm>
#include <cmath>

using namespace Rcpp;

//' @title Format a network table
//'
//' @param network_table Network edge table.
//' @param regulators,targets Nodes to include.
//' @param abs_weight Whether to use absolute weights and add interaction signs.
//' @return A formatted network edge table.
//' @export
// [[Rcpp::export]]
DataFrame network_format(DataFrame network_table,
                         Nullable<CharacterVector> regulators = R_NilValue,
                         Nullable<CharacterVector> targets = R_NilValue,
                         bool abs_weight = true)
{
  CharacterVector regulator = network_table["regulator"];
  CharacterVector target = network_table["target"];
  NumericVector weight = network_table["weight"];

  CharacterVector reg = regulators.isNotNull() ? CharacterVector(regulators) : CharacterVector();
  CharacterVector targ = targets.isNotNull() ? CharacterVector(targets) : CharacterVector();
  LogicalVector keep(weight.size(), false);
  for (int i = 0; i < weight.size(); ++i)
  {
    keep[i] = !NumericVector::is_na(weight[i]) && weight[i] != 0 &&
      (regulators.isNull() || std::find(reg.begin(), reg.end(), regulator[i]) != reg.end()) &&
      (targets.isNull() || std::find(targ.begin(), targ.end(), target[i]) != targ.end());
  }
  regulator = regulator[keep];
  target = target[keep];
  weight = weight[keep];

  CharacterVector interaction;
  if (abs_weight)
  {
    interaction = CharacterVector(weight.size());
    for (int i = 0; i < weight.size(); i++)
    {
      if (weight[i] < 0)
      {
        interaction[i] = "Repression";
        weight[i] = std::abs(weight[i]);
      }
      else
      {
        interaction[i] = "Activation";
      }
    }
  }

  IntegerVector order(weight.size());
  for (int i = 0; i < order.size(); i++)
    order[i] = i;
  std::sort(order.begin(), order.end(), [&](int i, int j) {
    return std::abs(weight[i]) > std::abs(weight[j]);
  });

  regulator = regulator[order];
  target = target[order];
  weight = weight[order];
  if (abs_weight)
  {
    interaction = interaction[order];
  }

  List result = List::create(
      Rcpp::Named("regulator") = regulator,
      Rcpp::Named("target") = target,
      Rcpp::Named("weight") = weight);
  if (abs_weight)
  {
    result.push_back(interaction, "Interaction");
  }

  return DataFrame(result);
}
