#include <Rcpp.h>
#include <unordered_map>
#include <string>

using namespace Rcpp;

// [[Rcpp::export]]
DataFrame prepare_calculate_metrics(DataFrame network_table, DataFrame ground_truth)
{
    CharacterVector network_reg = as<CharacterVector>(network_table["regulator"]);
    CharacterVector network_tar = as<CharacterVector>(network_table["target"]);
    NumericVector network_weight = abs(as<NumericVector>(network_table["weight"]));

    CharacterVector truth_reg = as<CharacterVector>(ground_truth["regulator"]);
    CharacterVector truth_tar = as<CharacterVector>(ground_truth["target"]);

    std::unordered_map<std::string, bool> truth_map;
    for (int i = 0; i < truth_reg.length(); i++)
    {
        if (truth_reg[i] == NA_STRING || truth_tar[i] == NA_STRING)
            continue;

        std::string key = as<std::string>(truth_reg[i]) + "|||" +
                          as<std::string>(truth_tar[i]);
        truth_map[key] = true;
    }

    int n = network_reg.length();
    CharacterVector out_reg(n);
    CharacterVector out_tar(n);
    NumericVector out_weight(n);
    IntegerVector out_label(n);

    for (int i = 0; i < n; i++)
    {
        if (network_reg[i] == NA_STRING || network_tar[i] == NA_STRING)
        {
            out_reg[i] = network_reg[i];
            out_tar[i] = network_tar[i];
            out_weight[i] = network_weight[i];
            out_label[i] = 0;
            continue;
        }

        std::string key = as<std::string>(network_reg[i]) + "|||" +
                          as<std::string>(network_tar[i]);

        out_reg[i] = network_reg[i];
        out_tar[i] = network_tar[i];
        out_weight[i] = network_weight[i];
        out_label[i] = truth_map.count(key) > 0 ? 1 : 0;
    }

    DataFrame result = DataFrame::create(
        _["regulator"] = out_reg,
        _["target"] = out_tar,
        _["weight"] = out_weight,
        _["label"] = out_label);

    IntegerVector row_names = seq(1, n);
    result.attr("row.names") = row_names;

    return result;
}

/*
prepare_calculate_metrics <- function(
    network_table,
    ground_truth) {
  colnames(network_table) <- c("regulator", "target", "weight")
  network_table$weight <- abs(as.numeric(network_table$weight))

  if (ncol(ground_truth) > 2) {
    ground_truth <- ground_truth[, 1:2]
  }
  names(ground_truth) <- c("regulator", "target")
  ground_truth$label <- rep(1, nrow(ground_truth))

  gold <- suppressWarnings(
    merge(
      network_table,
      ground_truth,
      by = c("regulator", "target"),
      all.x = TRUE
    )
  )
  gold$label[is.na(gold$label)] <- 0
  gold <- gold[order(gold$weight, decreasing = TRUE), ]
  rownames(gold) <- NULL

  return(gold)
}
*/