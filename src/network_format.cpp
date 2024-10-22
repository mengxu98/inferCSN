#include <Rcpp.h>
#include <algorithm>
#include <cmath>

using namespace Rcpp;

struct AbsGreater
{
  const NumericVector &weight;
  AbsGreater(const NumericVector &w) : weight(w) {}
  bool operator()(int i, int j) const
  {
    return std::abs(weight[i]) > std::abs(weight[j]);
  }
};

//' @title Format network table
//'
//' @param network_table The weight data table of network.
//' @param regulators Regulators list.
//' @param targets Targets list.
//' @param abs_weight Logical value, default is *`TRUE`*,
//' whether to perform absolute value on weights,
//' and when set `abs_weight` to *`TRUE`*,
//' the output of weight table will create a new column named `Interaction`.
//'
//' @md
//' @return Formated network table
//' @export
//'
//' @examples
//' data("example_matrix")
//' network_table <- inferCSN(example_matrix)
//'
//' network_format(
//'   network_table,
//'   regulators = "g1"
//' )
//'
//' network_format(
//'   network_table,
//'   regulators = "g1",
//'   abs_weight = FALSE
//' )
//'
//' network_format(
//'   network_table,
//'   targets = "g3"
//' )
//'
//' network_format(
//'   network_table,
//'   regulators = c("g1", "g3"),
//'   targets = c("g3", "g5")
//' )
// [[Rcpp::export]]
DataFrame network_format(DataFrame network_table,
                         Nullable<CharacterVector> regulators = R_NilValue,
                         Nullable<CharacterVector> targets = R_NilValue,
                         bool abs_weight = true)
{
  CharacterVector regulator = network_table["regulator"];
  CharacterVector target = network_table["target"];
  NumericVector weight = network_table["weight"];

  // filter rows with weight not equal to 0
  LogicalVector non_zero = weight != 0;
  regulator = regulator[non_zero];
  target = target[non_zero];
  weight = weight[non_zero];

  // process regulators
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

  // process targets
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

  // process abs_weight
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

  // sort
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

  // create result DataFrame
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

/*
#' @title Format network table
#'
#' @param network_table The weight data table of network.
#' @param regulators Regulators list.
#' @param targets Targets list.
#' @param abs_weight Logical value, default is *`TRUE`*,
#' whether to perform absolute value on weights,
#' and when set `abs_weight` to *`TRUE`*,
#' the output of weight table will create a new column named `Interaction`.
#'
#' @md
#' @return Formated network table
#' @export
#'
#' @examples
#' data("example_matrix")
#' network_table <- inferCSN(example_matrix)
#'
#' network_format(
#'   network_table,
#'   regulators = "g1"
#' )
#'
#' network_format(
#'   network_table,
#'   regulators = "g1",
#'   abs_weight = FALSE
#' )
#'
#' network_format(
#'   network_table,
#'   targets = "g3"
#' )
#'
#' network_format(
#'   network_table,
#'   regulators = c("g1", "g3"),
#'   targets = c("g3", "g5")
#' )
network_format <- function(
    network_table,
    regulators = NULL,
    targets = NULL,
    abs_weight = TRUE) {
  colnames(network_table) <- c("regulator", "target", "weight")
  network_table$weight <- as.numeric(network_table$weight)
  network_table <- dplyr::filter(network_table, weight != 0)
  if (!is.null(regulators)) {
    network_table <- purrr::map_dfr(
      unique(regulators), function(x) {
        dplyr::filter(network_table, regulator == x)
      }
    )
  }
  if (!is.null(targets)) {
    network_table <- purrr::map_dfr(
      unique(targets), function(x) {
        dplyr::filter(network_table, target == x)
      }
    )
  }

  if (abs_weight) {
    network_table$Interaction <- ifelse(
      network_table$weight < 0, "Repression", "Activation"
    )
    network_table$weight <- abs(network_table$weight)
  }

  network_table <- network_table[order(
    abs(as.numeric(network_table$weight)),
    decreasing = TRUE
  ), ]
  rownames(network_table) <- NULL

  return(network_table)
}
*/