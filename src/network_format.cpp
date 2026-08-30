#include "network_format.h"
#include <algorithm>
#include <cmath>

using namespace Rcpp;

AbsGreater::AbsGreater(const NumericVector &w) : weight(w) {}

bool AbsGreater::operator()(int i, int j) const
{
  return std::abs(weight[i]) > std::abs(weight[j]);
}

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

  LogicalVector not_na = !is_na(weight);
  regulator = regulator[not_na];
  target = target[not_na];
  weight = weight[not_na];

  LogicalVector non_zero = weight != 0;
  regulator = regulator[non_zero];
  target = target[non_zero];
  weight = weight[non_zero];

  if (regulators.isNotNull())
  {
    CharacterVector reg(regulators);
    LogicalVector keep(regulator.size(), false);
    for (int i = 0; i < regulator.size(); i++)
    {
      if (std::find(reg.begin(), reg.end(), regulator[i]) != reg.end())
      {
        keep[i] = true;
      }
    }
    regulator = regulator[keep];
    target = target[keep];
    weight = weight[keep];
  }

  if (targets.isNotNull())
  {
    CharacterVector targ(targets);
    LogicalVector keep(target.size(), false);
    for (int i = 0; i < target.size(); i++)
    {
      if (std::find(targ.begin(), targ.end(), target[i]) != targ.end())
      {
        keep[i] = true;
      }
    }
    regulator = regulator[keep];
    target = target[keep];
    weight = weight[keep];
  }

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
  std::sort(order.begin(), order.end(), AbsGreater(weight));

  regulator = regulator[order];
  target = target[order];
  weight = weight[order];
  if (abs_weight)
  {
    interaction = interaction[order];
  }

  DataFrame result;
  if (abs_weight)
  {
    result = DataFrame::create(
        Rcpp::Named("regulator") = regulator,
        Rcpp::Named("target") = target,
        Rcpp::Named("weight") = weight,
        Rcpp::Named("Interaction") = interaction);
  }
  else
  {
    result = DataFrame::create(
        Rcpp::Named("regulator") = regulator,
        Rcpp::Named("target") = target,
        Rcpp::Named("weight") = weight);
  }

  return result;
}
