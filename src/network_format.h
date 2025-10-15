#ifndef NETWORK_FORMAT_H
#define NETWORK_FORMAT_H

#pragma once
#include <Rcpp.h>
using namespace Rcpp;

struct AbsGreater
{
  const Rcpp::NumericVector &weight;
  AbsGreater(const Rcpp::NumericVector &w);
  bool operator()(int i, int j) const;
};

Rcpp::DataFrame network_format(
    Rcpp::DataFrame network_table,
    Rcpp::Nullable<Rcpp::CharacterVector> regulators,
    Rcpp::Nullable<Rcpp::CharacterVector> targets,
    bool abs_weight);

#endif
