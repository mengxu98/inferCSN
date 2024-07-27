#ifndef PARAMS_H
#define PARAMS_H
#include "BetaVector.h"
#include "Model.h"
#include "RcppArmadillo.h"
#include <map>

// Version 1
// template <typename T>
// struct Params {

//     Model Specs;
//     std::vector<double> ModelParams {0, 0, 0, 2};
//     std::size_t MaxIters = 500;
//     double rtol = 1e-8;
//     double atol = 1e-12;
//     char Init = 'z'; // 'z' => zeros
//     std::size_t RandomStartSize = 10;

//     // beta_vector * InitialSol; // This line code will encounter warnings
//     following:
//       // Found the following significant warnings:
//         // include/Params.h:9:8: warning: 'P.Params<arma::SpMat<double>
//         >::InitialSol' is used uninitialized [-Wuninitialized]
//         // include/Params.h:9:8: warning: 'P' is used uninitialized
//         [-Wuninitialized]
//         // include/Params.h:9:8: warning: 'P.Params<arma::Mat<double>
//         >::InitialSol' is used uninitialized [-Wuninitialized]
//     // First, I use code: 'beta_vector * InitialSol = nullptr' to attempt to
//     eliminate warnings
//     // beta_vector * InitialSol = nullptr; // Initialize to nullptr,
//     indicating no initial solution
//     // But it also prompts:
//       // Found the following significant warnings:
//         // include/Params.h:9:8: warning: ‘P’ is used uninitialized
//         [-Wuninitialized]
//     // Now, I use code: 'beta_vector * InitialSol' to attempt again
//     arma::vec* InitialSol;

//     double b0 = 0; // intercept
//     char CyclingOrder = 'c';
//     std::vector<std::size_t> Uorder;
//     bool ActiveSet = true;
//     std::size_t ActiveSetNum = 6;
//     std::size_t MaxNumSwaps = 200; // Used by CDSwaps
//     std::vector<double> * Xtr;
//     arma::rowvec * ytX;
//     std::map<std::size_t, arma::rowvec> * D;
//     std::size_t Iter = 0; // Current iteration number in the grid
//     std::size_t ScreenSize = 1000;
//     arma::vec * r;
//     T * Xy; // used for classification.
//     std::size_t NoSelectK = 0;
//     bool intercept = false;
//     bool withBounds = false;
//     arma::vec Lows;
//     arma::vec Highs;

// };

// #endif

// Version 2
template <typename T> struct Params {

  Model Specs;
  std::vector<double> ModelParams{0, 0, 0, 2};
  std::size_t MaxIters = 500;
  double rtol = 1e-8;
  double atol = 1e-12;
  char Init = 'z';
  std::size_t RandomStartSize = 10;

  // beta_vector * InitialSol; // This line code will encounter warnings
  // following:
  //   // Found the following significant warnings:
  //     // include/Params.h:9:8: warning: 'P.Params<arma::SpMat<double>
  //     >::InitialSol' is used uninitialized [-Wuninitialized]
  //     // include/Params.h:9:8: warning: 'P' is used uninitialized
  //     [-Wuninitialized]
  //     // include/Params.h:9:8: warning: 'P.Params<arma::Mat<double>
  //     >::InitialSol' is used uninitialized [-Wuninitialized]
  // // First, I use code: 'beta_vector * InitialSol = nullptr' to attempt to
  // eliminate warnings beta_vector * InitialSol = nullptr; // Initialize to
  // nullptr, indicating no initial solution
  //     // But it also prompts:
  //         // Found the following significant warnings:
  //         // include/Params.h:9:8: warning: ‘P’ is used uninitialized
  //         [-Wuninitialized]
  //         // Now, I use code: 'beta_vector * InitialSol' to attempt again
  // arma::vec* InitialSol;
  arma::vec *InitialSol = nullptr;

  double b0 = 0;
  char CyclingOrder = 'c';
  std::vector<std::size_t> Uorder;
  bool ActiveSet = true;
  std::size_t ActiveSetNum = 6;
  std::size_t MaxNumSwaps = 200;
  std::vector<double> *Xtr = nullptr;
  arma::rowvec *ytX = nullptr;
  std::map<std::size_t, arma::rowvec> *D = nullptr;
  std::size_t Iter = 0;
  std::size_t ScreenSize = 1000;
  arma::vec *r = nullptr;
  T *Xy = nullptr;
  std::size_t NoSelectK = 0;
  bool intercept = false;
  bool withBounds = false;
  arma::vec Lows;
  arma::vec Highs;

  // Params() : InitialSol(nullptr) { }
  Params(){};
};

#endif
