#ifndef GRIDPARAMS_H
#define GRIDPARAMS_H
#include "Params.h"
#include "RcppArmadillo.h"

template <typename T> struct GridParams {

  Params<T> P;
  std::size_t G_ncols = 100;
  std::size_t G_nrows = 10;
  bool LambdaU = false;
  std::size_t NnzStopNum = 200;
  double LambdaMinFactor = 0.01;
  arma::vec Lambdas;
  std::vector<std::vector<double>> LambdasGrid;
  double Lambda2Max = 0.1;
  double Lambda2Min = 0.001;
  std::string Type = "L0";
  bool PartialSort = true;
  bool XtrAvailable = false;
  double ytXmax;
  std::vector<double> *Xtr;
  double ScaleDownFactor = 0.8;
  bool intercept;

  GridParams(){}; // Used to fix errors when compiling C++source code, as
                  // follows:
                  //     include/GridParams.h: In instantiation of
                  //     ‘Grid<T>::Grid(const T&, const vec&, const
                  //     GridParams<T>&) [with T = arma::Mat<double>; arma::vec
                  //     = arma::Col<double>]’:
                  // Grid.cpp:69:16:   required from here
                  // include/GridParams.h:7:8: error: conversion from ‘const
                  // char [3]’ to non-scalar type ‘std::string {aka
                  // std::basic_string<char>}’ requested
                  //  struct GridParams
};

#endif
