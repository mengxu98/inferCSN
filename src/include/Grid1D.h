#ifndef GRID1D_H
#define GRID1D_H
#include "FitResult.h"
#include "GridParams.h"
#include "MakeCD.h"
#include "Params.h"
#include "RcppArmadillo.h"
#include <algorithm>
#include <map>
#include <memory>

template <class T> class Grid1D {
private:
  std::size_t G_ncols;
  Params<T> P;
  const T *X;
  const arma::vec *y;
  std::size_t p;
  std::vector<std::unique_ptr<FitResult<T>>> G;
  arma::vec Lambdas;
  bool LambdaU;
  std::size_t NnzStopNum;
  std::vector<double> *Xtr;
  arma::rowvec *ytX;
  double LambdaMinFactor;
  bool PartialSort;
  bool XtrAvailable;
  double ytXmax2d;
  double ScaleDownFactor;
  std::size_t NoSelectK;

public:
  Grid1D(const T &Xi, const arma::vec &yi, const GridParams<T> &PG);
  ~Grid1D();
  std::vector<std::unique_ptr<FitResult<T>>> Fit();
};

#endif
