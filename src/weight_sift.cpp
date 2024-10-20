#include <Rcpp.h>
#include <unordered_set>
#include <algorithm>
#include <cmath>

//' @title Weight sift
//' @description Remove edges with smaller weights in the reverse direction.
//'
//' @param table A data frame with three columns: "regulator", "target", and "weight".
//'
//' @export
//'
//' @examples
//' data("example_matrix")
//' network_table <- inferCSN(example_matrix)
//' weight_sift(network_table) |> head()
// [[Rcpp::export]]
Rcpp::DataFrame weight_sift(Rcpp::DataFrame table) {
    Rcpp::StringVector x = table["regulator"];
    Rcpp::StringVector y = table["target"];
    Rcpp::NumericVector v = table["weight"];
    
    int n = x.size();
    
    // Preallocate memory
    Rcpp::StringVector edges(n);
    Rcpp::StringVector reversed_edges(n);
    Rcpp::NumericVector reversed_v(n);
    
    std::string edge, reversed_edge;
    edge.reserve(50);  // Estimated maximum length
    reversed_edge.reserve(50);
    
    for(int i = 0; i < n; ++i) {
        edge = x[i];
        edge += "_";
        edge += y[i];
        edges[i] = edge;
        
        reversed_edge = y[i];
        reversed_edge += "_";
        reversed_edge += x[i];
        reversed_edges[i] = reversed_edge;
        reversed_v[i] = v[i];
    }
    
    // Use Rcpp sugar functions
    Rcpp::LogicalVector is_common = Rcpp::in(edges, reversed_edges);
    
    if(Rcpp::sum(is_common) == 0) return table;
    
    Rcpp::LogicalVector keep(n, true);
    Rcpp::IntegerVector match_indices = Rcpp::match(edges, reversed_edges) - 1;
    
    for(int i = 0; i < n; ++i) {
        if(is_common[i]) {
            int rev_index = match_indices[i];
            if(std::abs(v[i]) < std::abs(reversed_v[rev_index])) {
                keep[i] = false;
            }
        }
    }
    
    return Rcpp::DataFrame::create(
        Rcpp::_["regulator"] = x[keep],
        Rcpp::_["target"] = y[keep],
        Rcpp::_["weight"] = v[keep]
    );
}

/*
# R code
weight_sift <- function(table) {
  table <- table[, 1:3]
  raw_rownames <- colnames(table)
  colnames(table) <- c("x", "y", "v")
  table$edge <- paste(
    table$x,
    table$y,
    sep = "_"
  )
  rownames(table) <- table$edge

  table_new <- data.frame(
    edge = paste(
      table$y,
      table$x,
      sep = "_"
    ),
    v = table$v
  )
  rownames(table_new) <- table_new$edge

  common_edges <- intersect(rownames(table), rownames(table_new))
  if (length(common_edges) == 0) {
    table <- table[, 1:3]
    colnames(table) <- raw_rownames
    rownames(table) <- NULL
    return(table)
  }
  table_common_edges <- table[common_edges, ]
  table_new <- table_new[common_edges, ]
  table_common_edges$v_new <- table_new$v
  table_common_edges <- dplyr::filter(table_common_edges, abs(v) < abs(v_new))
  table <- table[setdiff(rownames(table), rownames(table_common_edges)), 1:3]
  colnames(table) <- raw_rownames
  rownames(table) <- NULL

  return(table)
}
*/