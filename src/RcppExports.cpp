// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// srm_model_sparse
Rcpp::List srm_model_sparse(const arma::sp_mat& X, const arma::vec& y, const std::string Loss, const std::string Penalty, const std::string Algorithm, const std::size_t NnzStopNum, const std::size_t G_ncols, const std::size_t G_nrows, const double Lambda2Max, const double Lambda2Min, const bool PartialSort, const std::size_t MaxIters, const double rtol, const double atol, const bool ActiveSet, const std::size_t ActiveSetNum, const std::size_t MaxNumSwaps, const double ScaleDownFactor, const std::size_t ScreenSize, const bool LambdaU, const std::vector<std::vector<double>> Lambdas, const std::size_t ExcludeFirstK, const bool Intercept, const bool withBounds, const arma::vec& Lows, const arma::vec& Highs);
RcppExport SEXP _inferCSN_srm_model_sparse(SEXP XSEXP, SEXP ySEXP, SEXP LossSEXP, SEXP PenaltySEXP, SEXP AlgorithmSEXP, SEXP NnzStopNumSEXP, SEXP G_ncolsSEXP, SEXP G_nrowsSEXP, SEXP Lambda2MaxSEXP, SEXP Lambda2MinSEXP, SEXP PartialSortSEXP, SEXP MaxItersSEXP, SEXP rtolSEXP, SEXP atolSEXP, SEXP ActiveSetSEXP, SEXP ActiveSetNumSEXP, SEXP MaxNumSwapsSEXP, SEXP ScaleDownFactorSEXP, SEXP ScreenSizeSEXP, SEXP LambdaUSEXP, SEXP LambdasSEXP, SEXP ExcludeFirstKSEXP, SEXP InterceptSEXP, SEXP withBoundsSEXP, SEXP LowsSEXP, SEXP HighsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const std::string >::type Loss(LossSEXP);
    Rcpp::traits::input_parameter< const std::string >::type Penalty(PenaltySEXP);
    Rcpp::traits::input_parameter< const std::string >::type Algorithm(AlgorithmSEXP);
    Rcpp::traits::input_parameter< const std::size_t >::type NnzStopNum(NnzStopNumSEXP);
    Rcpp::traits::input_parameter< const std::size_t >::type G_ncols(G_ncolsSEXP);
    Rcpp::traits::input_parameter< const std::size_t >::type G_nrows(G_nrowsSEXP);
    Rcpp::traits::input_parameter< const double >::type Lambda2Max(Lambda2MaxSEXP);
    Rcpp::traits::input_parameter< const double >::type Lambda2Min(Lambda2MinSEXP);
    Rcpp::traits::input_parameter< const bool >::type PartialSort(PartialSortSEXP);
    Rcpp::traits::input_parameter< const std::size_t >::type MaxIters(MaxItersSEXP);
    Rcpp::traits::input_parameter< const double >::type rtol(rtolSEXP);
    Rcpp::traits::input_parameter< const double >::type atol(atolSEXP);
    Rcpp::traits::input_parameter< const bool >::type ActiveSet(ActiveSetSEXP);
    Rcpp::traits::input_parameter< const std::size_t >::type ActiveSetNum(ActiveSetNumSEXP);
    Rcpp::traits::input_parameter< const std::size_t >::type MaxNumSwaps(MaxNumSwapsSEXP);
    Rcpp::traits::input_parameter< const double >::type ScaleDownFactor(ScaleDownFactorSEXP);
    Rcpp::traits::input_parameter< const std::size_t >::type ScreenSize(ScreenSizeSEXP);
    Rcpp::traits::input_parameter< const bool >::type LambdaU(LambdaUSEXP);
    Rcpp::traits::input_parameter< const std::vector<std::vector<double>> >::type Lambdas(LambdasSEXP);
    Rcpp::traits::input_parameter< const std::size_t >::type ExcludeFirstK(ExcludeFirstKSEXP);
    Rcpp::traits::input_parameter< const bool >::type Intercept(InterceptSEXP);
    Rcpp::traits::input_parameter< const bool >::type withBounds(withBoundsSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type Lows(LowsSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type Highs(HighsSEXP);
    rcpp_result_gen = Rcpp::wrap(srm_model_sparse(X, y, Loss, Penalty, Algorithm, NnzStopNum, G_ncols, G_nrows, Lambda2Max, Lambda2Min, PartialSort, MaxIters, rtol, atol, ActiveSet, ActiveSetNum, MaxNumSwaps, ScaleDownFactor, ScreenSize, LambdaU, Lambdas, ExcludeFirstK, Intercept, withBounds, Lows, Highs));
    return rcpp_result_gen;
END_RCPP
}
// srm_model_dense
Rcpp::List srm_model_dense(const arma::mat& X, const arma::vec& y, const std::string Loss, const std::string Penalty, const std::string Algorithm, const std::size_t NnzStopNum, const std::size_t G_ncols, const std::size_t G_nrows, const double Lambda2Max, const double Lambda2Min, const bool PartialSort, const std::size_t MaxIters, const double rtol, const double atol, const bool ActiveSet, const std::size_t ActiveSetNum, const std::size_t MaxNumSwaps, const double ScaleDownFactor, const std::size_t ScreenSize, const bool LambdaU, const std::vector<std::vector<double>> Lambdas, const std::size_t ExcludeFirstK, const bool Intercept, const bool withBounds, const arma::vec& Lows, const arma::vec& Highs);
RcppExport SEXP _inferCSN_srm_model_dense(SEXP XSEXP, SEXP ySEXP, SEXP LossSEXP, SEXP PenaltySEXP, SEXP AlgorithmSEXP, SEXP NnzStopNumSEXP, SEXP G_ncolsSEXP, SEXP G_nrowsSEXP, SEXP Lambda2MaxSEXP, SEXP Lambda2MinSEXP, SEXP PartialSortSEXP, SEXP MaxItersSEXP, SEXP rtolSEXP, SEXP atolSEXP, SEXP ActiveSetSEXP, SEXP ActiveSetNumSEXP, SEXP MaxNumSwapsSEXP, SEXP ScaleDownFactorSEXP, SEXP ScreenSizeSEXP, SEXP LambdaUSEXP, SEXP LambdasSEXP, SEXP ExcludeFirstKSEXP, SEXP InterceptSEXP, SEXP withBoundsSEXP, SEXP LowsSEXP, SEXP HighsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const std::string >::type Loss(LossSEXP);
    Rcpp::traits::input_parameter< const std::string >::type Penalty(PenaltySEXP);
    Rcpp::traits::input_parameter< const std::string >::type Algorithm(AlgorithmSEXP);
    Rcpp::traits::input_parameter< const std::size_t >::type NnzStopNum(NnzStopNumSEXP);
    Rcpp::traits::input_parameter< const std::size_t >::type G_ncols(G_ncolsSEXP);
    Rcpp::traits::input_parameter< const std::size_t >::type G_nrows(G_nrowsSEXP);
    Rcpp::traits::input_parameter< const double >::type Lambda2Max(Lambda2MaxSEXP);
    Rcpp::traits::input_parameter< const double >::type Lambda2Min(Lambda2MinSEXP);
    Rcpp::traits::input_parameter< const bool >::type PartialSort(PartialSortSEXP);
    Rcpp::traits::input_parameter< const std::size_t >::type MaxIters(MaxItersSEXP);
    Rcpp::traits::input_parameter< const double >::type rtol(rtolSEXP);
    Rcpp::traits::input_parameter< const double >::type atol(atolSEXP);
    Rcpp::traits::input_parameter< const bool >::type ActiveSet(ActiveSetSEXP);
    Rcpp::traits::input_parameter< const std::size_t >::type ActiveSetNum(ActiveSetNumSEXP);
    Rcpp::traits::input_parameter< const std::size_t >::type MaxNumSwaps(MaxNumSwapsSEXP);
    Rcpp::traits::input_parameter< const double >::type ScaleDownFactor(ScaleDownFactorSEXP);
    Rcpp::traits::input_parameter< const std::size_t >::type ScreenSize(ScreenSizeSEXP);
    Rcpp::traits::input_parameter< const bool >::type LambdaU(LambdaUSEXP);
    Rcpp::traits::input_parameter< const std::vector<std::vector<double>> >::type Lambdas(LambdasSEXP);
    Rcpp::traits::input_parameter< const std::size_t >::type ExcludeFirstK(ExcludeFirstKSEXP);
    Rcpp::traits::input_parameter< const bool >::type Intercept(InterceptSEXP);
    Rcpp::traits::input_parameter< const bool >::type withBounds(withBoundsSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type Lows(LowsSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type Highs(HighsSEXP);
    rcpp_result_gen = Rcpp::wrap(srm_model_dense(X, y, Loss, Penalty, Algorithm, NnzStopNum, G_ncols, G_nrows, Lambda2Max, Lambda2Min, PartialSort, MaxIters, rtol, atol, ActiveSet, ActiveSetNum, MaxNumSwaps, ScaleDownFactor, ScreenSize, LambdaU, Lambdas, ExcludeFirstK, Intercept, withBounds, Lows, Highs));
    return rcpp_result_gen;
END_RCPP
}
// srm_model_cv_sparse
Rcpp::List srm_model_cv_sparse(const arma::sp_mat& X, const arma::vec& y, const std::string Loss, const std::string Penalty, const std::string Algorithm, const std::size_t NnzStopNum, const std::size_t G_ncols, const std::size_t G_nrows, const double Lambda2Max, const double Lambda2Min, const bool PartialSort, const std::size_t MaxIters, const double rtol, const double atol, const bool ActiveSet, const std::size_t ActiveSetNum, const std::size_t MaxNumSwaps, const double ScaleDownFactor, const std::size_t ScreenSize, const bool LambdaU, const std::vector<std::vector<double>> Lambdas, const std::size_t nfolds, const double seed, const std::size_t ExcludeFirstK, const bool Intercept, const bool withBounds, const arma::vec& Lows, const arma::vec& Highs);
RcppExport SEXP _inferCSN_srm_model_cv_sparse(SEXP XSEXP, SEXP ySEXP, SEXP LossSEXP, SEXP PenaltySEXP, SEXP AlgorithmSEXP, SEXP NnzStopNumSEXP, SEXP G_ncolsSEXP, SEXP G_nrowsSEXP, SEXP Lambda2MaxSEXP, SEXP Lambda2MinSEXP, SEXP PartialSortSEXP, SEXP MaxItersSEXP, SEXP rtolSEXP, SEXP atolSEXP, SEXP ActiveSetSEXP, SEXP ActiveSetNumSEXP, SEXP MaxNumSwapsSEXP, SEXP ScaleDownFactorSEXP, SEXP ScreenSizeSEXP, SEXP LambdaUSEXP, SEXP LambdasSEXP, SEXP nfoldsSEXP, SEXP seedSEXP, SEXP ExcludeFirstKSEXP, SEXP InterceptSEXP, SEXP withBoundsSEXP, SEXP LowsSEXP, SEXP HighsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const std::string >::type Loss(LossSEXP);
    Rcpp::traits::input_parameter< const std::string >::type Penalty(PenaltySEXP);
    Rcpp::traits::input_parameter< const std::string >::type Algorithm(AlgorithmSEXP);
    Rcpp::traits::input_parameter< const std::size_t >::type NnzStopNum(NnzStopNumSEXP);
    Rcpp::traits::input_parameter< const std::size_t >::type G_ncols(G_ncolsSEXP);
    Rcpp::traits::input_parameter< const std::size_t >::type G_nrows(G_nrowsSEXP);
    Rcpp::traits::input_parameter< const double >::type Lambda2Max(Lambda2MaxSEXP);
    Rcpp::traits::input_parameter< const double >::type Lambda2Min(Lambda2MinSEXP);
    Rcpp::traits::input_parameter< const bool >::type PartialSort(PartialSortSEXP);
    Rcpp::traits::input_parameter< const std::size_t >::type MaxIters(MaxItersSEXP);
    Rcpp::traits::input_parameter< const double >::type rtol(rtolSEXP);
    Rcpp::traits::input_parameter< const double >::type atol(atolSEXP);
    Rcpp::traits::input_parameter< const bool >::type ActiveSet(ActiveSetSEXP);
    Rcpp::traits::input_parameter< const std::size_t >::type ActiveSetNum(ActiveSetNumSEXP);
    Rcpp::traits::input_parameter< const std::size_t >::type MaxNumSwaps(MaxNumSwapsSEXP);
    Rcpp::traits::input_parameter< const double >::type ScaleDownFactor(ScaleDownFactorSEXP);
    Rcpp::traits::input_parameter< const std::size_t >::type ScreenSize(ScreenSizeSEXP);
    Rcpp::traits::input_parameter< const bool >::type LambdaU(LambdaUSEXP);
    Rcpp::traits::input_parameter< const std::vector<std::vector<double>> >::type Lambdas(LambdasSEXP);
    Rcpp::traits::input_parameter< const std::size_t >::type nfolds(nfoldsSEXP);
    Rcpp::traits::input_parameter< const double >::type seed(seedSEXP);
    Rcpp::traits::input_parameter< const std::size_t >::type ExcludeFirstK(ExcludeFirstKSEXP);
    Rcpp::traits::input_parameter< const bool >::type Intercept(InterceptSEXP);
    Rcpp::traits::input_parameter< const bool >::type withBounds(withBoundsSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type Lows(LowsSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type Highs(HighsSEXP);
    rcpp_result_gen = Rcpp::wrap(srm_model_cv_sparse(X, y, Loss, Penalty, Algorithm, NnzStopNum, G_ncols, G_nrows, Lambda2Max, Lambda2Min, PartialSort, MaxIters, rtol, atol, ActiveSet, ActiveSetNum, MaxNumSwaps, ScaleDownFactor, ScreenSize, LambdaU, Lambdas, nfolds, seed, ExcludeFirstK, Intercept, withBounds, Lows, Highs));
    return rcpp_result_gen;
END_RCPP
}
// srm_model_cv_dense
Rcpp::List srm_model_cv_dense(const arma::mat& X, const arma::vec& y, const std::string Loss, const std::string Penalty, const std::string Algorithm, const std::size_t NnzStopNum, const std::size_t G_ncols, const std::size_t G_nrows, const double Lambda2Max, const double Lambda2Min, const bool PartialSort, const std::size_t MaxIters, const double rtol, const double atol, const bool ActiveSet, const std::size_t ActiveSetNum, const std::size_t MaxNumSwaps, const double ScaleDownFactor, const std::size_t ScreenSize, const bool LambdaU, const std::vector<std::vector<double>> Lambdas, const std::size_t nfolds, const double seed, const std::size_t ExcludeFirstK, const bool Intercept, const bool withBounds, const arma::vec& Lows, const arma::vec& Highs);
RcppExport SEXP _inferCSN_srm_model_cv_dense(SEXP XSEXP, SEXP ySEXP, SEXP LossSEXP, SEXP PenaltySEXP, SEXP AlgorithmSEXP, SEXP NnzStopNumSEXP, SEXP G_ncolsSEXP, SEXP G_nrowsSEXP, SEXP Lambda2MaxSEXP, SEXP Lambda2MinSEXP, SEXP PartialSortSEXP, SEXP MaxItersSEXP, SEXP rtolSEXP, SEXP atolSEXP, SEXP ActiveSetSEXP, SEXP ActiveSetNumSEXP, SEXP MaxNumSwapsSEXP, SEXP ScaleDownFactorSEXP, SEXP ScreenSizeSEXP, SEXP LambdaUSEXP, SEXP LambdasSEXP, SEXP nfoldsSEXP, SEXP seedSEXP, SEXP ExcludeFirstKSEXP, SEXP InterceptSEXP, SEXP withBoundsSEXP, SEXP LowsSEXP, SEXP HighsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const std::string >::type Loss(LossSEXP);
    Rcpp::traits::input_parameter< const std::string >::type Penalty(PenaltySEXP);
    Rcpp::traits::input_parameter< const std::string >::type Algorithm(AlgorithmSEXP);
    Rcpp::traits::input_parameter< const std::size_t >::type NnzStopNum(NnzStopNumSEXP);
    Rcpp::traits::input_parameter< const std::size_t >::type G_ncols(G_ncolsSEXP);
    Rcpp::traits::input_parameter< const std::size_t >::type G_nrows(G_nrowsSEXP);
    Rcpp::traits::input_parameter< const double >::type Lambda2Max(Lambda2MaxSEXP);
    Rcpp::traits::input_parameter< const double >::type Lambda2Min(Lambda2MinSEXP);
    Rcpp::traits::input_parameter< const bool >::type PartialSort(PartialSortSEXP);
    Rcpp::traits::input_parameter< const std::size_t >::type MaxIters(MaxItersSEXP);
    Rcpp::traits::input_parameter< const double >::type rtol(rtolSEXP);
    Rcpp::traits::input_parameter< const double >::type atol(atolSEXP);
    Rcpp::traits::input_parameter< const bool >::type ActiveSet(ActiveSetSEXP);
    Rcpp::traits::input_parameter< const std::size_t >::type ActiveSetNum(ActiveSetNumSEXP);
    Rcpp::traits::input_parameter< const std::size_t >::type MaxNumSwaps(MaxNumSwapsSEXP);
    Rcpp::traits::input_parameter< const double >::type ScaleDownFactor(ScaleDownFactorSEXP);
    Rcpp::traits::input_parameter< const std::size_t >::type ScreenSize(ScreenSizeSEXP);
    Rcpp::traits::input_parameter< const bool >::type LambdaU(LambdaUSEXP);
    Rcpp::traits::input_parameter< const std::vector<std::vector<double>> >::type Lambdas(LambdasSEXP);
    Rcpp::traits::input_parameter< const std::size_t >::type nfolds(nfoldsSEXP);
    Rcpp::traits::input_parameter< const double >::type seed(seedSEXP);
    Rcpp::traits::input_parameter< const std::size_t >::type ExcludeFirstK(ExcludeFirstKSEXP);
    Rcpp::traits::input_parameter< const bool >::type Intercept(InterceptSEXP);
    Rcpp::traits::input_parameter< const bool >::type withBounds(withBoundsSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type Lows(LowsSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type Highs(HighsSEXP);
    rcpp_result_gen = Rcpp::wrap(srm_model_cv_dense(X, y, Loss, Penalty, Algorithm, NnzStopNum, G_ncols, G_nrows, Lambda2Max, Lambda2Min, PartialSort, MaxIters, rtol, atol, ActiveSet, ActiveSetNum, MaxNumSwaps, ScaleDownFactor, ScreenSize, LambdaU, Lambdas, nfolds, seed, ExcludeFirstK, Intercept, withBounds, Lows, Highs));
    return rcpp_result_gen;
END_RCPP
}
// R_matrix_column_get_dense
arma::vec R_matrix_column_get_dense(const arma::mat& mat, int col);
RcppExport SEXP _inferCSN_R_matrix_column_get_dense(SEXP matSEXP, SEXP colSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type mat(matSEXP);
    Rcpp::traits::input_parameter< int >::type col(colSEXP);
    rcpp_result_gen = Rcpp::wrap(R_matrix_column_get_dense(mat, col));
    return rcpp_result_gen;
END_RCPP
}
// R_matrix_column_get_sparse
arma::vec R_matrix_column_get_sparse(const arma::sp_mat& mat, int col);
RcppExport SEXP _inferCSN_R_matrix_column_get_sparse(SEXP matSEXP, SEXP colSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type mat(matSEXP);
    Rcpp::traits::input_parameter< int >::type col(colSEXP);
    rcpp_result_gen = Rcpp::wrap(R_matrix_column_get_sparse(mat, col));
    return rcpp_result_gen;
END_RCPP
}
// R_matrix_rows_get_dense
arma::mat R_matrix_rows_get_dense(const arma::mat& mat, const arma::ucolvec rows);
RcppExport SEXP _inferCSN_R_matrix_rows_get_dense(SEXP matSEXP, SEXP rowsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type mat(matSEXP);
    Rcpp::traits::input_parameter< const arma::ucolvec >::type rows(rowsSEXP);
    rcpp_result_gen = Rcpp::wrap(R_matrix_rows_get_dense(mat, rows));
    return rcpp_result_gen;
END_RCPP
}
// R_matrix_rows_get_sparse
arma::sp_mat R_matrix_rows_get_sparse(const arma::sp_mat& mat, const arma::ucolvec rows);
RcppExport SEXP _inferCSN_R_matrix_rows_get_sparse(SEXP matSEXP, SEXP rowsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type mat(matSEXP);
    Rcpp::traits::input_parameter< const arma::ucolvec >::type rows(rowsSEXP);
    rcpp_result_gen = Rcpp::wrap(R_matrix_rows_get_sparse(mat, rows));
    return rcpp_result_gen;
END_RCPP
}
// R_matrix_vector_schur_product_dense
arma::mat R_matrix_vector_schur_product_dense(const arma::mat& mat, const arma::vec& u);
RcppExport SEXP _inferCSN_R_matrix_vector_schur_product_dense(SEXP matSEXP, SEXP uSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type mat(matSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type u(uSEXP);
    rcpp_result_gen = Rcpp::wrap(R_matrix_vector_schur_product_dense(mat, u));
    return rcpp_result_gen;
END_RCPP
}
// R_matrix_vector_schur_product_sparse
arma::sp_mat R_matrix_vector_schur_product_sparse(const arma::sp_mat& mat, const arma::vec& u);
RcppExport SEXP _inferCSN_R_matrix_vector_schur_product_sparse(SEXP matSEXP, SEXP uSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type mat(matSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type u(uSEXP);
    rcpp_result_gen = Rcpp::wrap(R_matrix_vector_schur_product_sparse(mat, u));
    return rcpp_result_gen;
END_RCPP
}
// R_matrix_vector_divide_dense
arma::mat R_matrix_vector_divide_dense(const arma::mat& mat, const arma::vec& u);
RcppExport SEXP _inferCSN_R_matrix_vector_divide_dense(SEXP matSEXP, SEXP uSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type mat(matSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type u(uSEXP);
    rcpp_result_gen = Rcpp::wrap(R_matrix_vector_divide_dense(mat, u));
    return rcpp_result_gen;
END_RCPP
}
// R_matrix_vector_divide_sparse
arma::sp_mat R_matrix_vector_divide_sparse(const arma::sp_mat& mat, const arma::vec& u);
RcppExport SEXP _inferCSN_R_matrix_vector_divide_sparse(SEXP matSEXP, SEXP uSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type mat(matSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type u(uSEXP);
    rcpp_result_gen = Rcpp::wrap(R_matrix_vector_divide_sparse(mat, u));
    return rcpp_result_gen;
END_RCPP
}
// R_matrix_column_sums_dense
arma::rowvec R_matrix_column_sums_dense(const arma::mat& mat);
RcppExport SEXP _inferCSN_R_matrix_column_sums_dense(SEXP matSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type mat(matSEXP);
    rcpp_result_gen = Rcpp::wrap(R_matrix_column_sums_dense(mat));
    return rcpp_result_gen;
END_RCPP
}
// R_matrix_column_sums_sparse
arma::rowvec R_matrix_column_sums_sparse(const arma::sp_mat& mat);
RcppExport SEXP _inferCSN_R_matrix_column_sums_sparse(SEXP matSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type mat(matSEXP);
    rcpp_result_gen = Rcpp::wrap(R_matrix_column_sums_sparse(mat));
    return rcpp_result_gen;
END_RCPP
}
// R_matrix_column_dot_dense
double R_matrix_column_dot_dense(const arma::mat& mat, int col, const arma::vec u);
RcppExport SEXP _inferCSN_R_matrix_column_dot_dense(SEXP matSEXP, SEXP colSEXP, SEXP uSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type mat(matSEXP);
    Rcpp::traits::input_parameter< int >::type col(colSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type u(uSEXP);
    rcpp_result_gen = Rcpp::wrap(R_matrix_column_dot_dense(mat, col, u));
    return rcpp_result_gen;
END_RCPP
}
// R_matrix_column_dot_sparse
double R_matrix_column_dot_sparse(const arma::sp_mat& mat, int col, const arma::vec u);
RcppExport SEXP _inferCSN_R_matrix_column_dot_sparse(SEXP matSEXP, SEXP colSEXP, SEXP uSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type mat(matSEXP);
    Rcpp::traits::input_parameter< int >::type col(colSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type u(uSEXP);
    rcpp_result_gen = Rcpp::wrap(R_matrix_column_dot_sparse(mat, col, u));
    return rcpp_result_gen;
END_RCPP
}
// R_matrix_column_mult_dense
arma::vec R_matrix_column_mult_dense(const arma::mat& mat, int col, double u);
RcppExport SEXP _inferCSN_R_matrix_column_mult_dense(SEXP matSEXP, SEXP colSEXP, SEXP uSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type mat(matSEXP);
    Rcpp::traits::input_parameter< int >::type col(colSEXP);
    Rcpp::traits::input_parameter< double >::type u(uSEXP);
    rcpp_result_gen = Rcpp::wrap(R_matrix_column_mult_dense(mat, col, u));
    return rcpp_result_gen;
END_RCPP
}
// R_matrix_column_mult_sparse
arma::vec R_matrix_column_mult_sparse(const arma::sp_mat& mat, int col, double u);
RcppExport SEXP _inferCSN_R_matrix_column_mult_sparse(SEXP matSEXP, SEXP colSEXP, SEXP uSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type mat(matSEXP);
    Rcpp::traits::input_parameter< int >::type col(colSEXP);
    Rcpp::traits::input_parameter< double >::type u(uSEXP);
    rcpp_result_gen = Rcpp::wrap(R_matrix_column_mult_sparse(mat, col, u));
    return rcpp_result_gen;
END_RCPP
}
// R_matrix_normalize_dense
Rcpp::List R_matrix_normalize_dense(arma::mat mat_norm);
RcppExport SEXP _inferCSN_R_matrix_normalize_dense(SEXP mat_normSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type mat_norm(mat_normSEXP);
    rcpp_result_gen = Rcpp::wrap(R_matrix_normalize_dense(mat_norm));
    return rcpp_result_gen;
END_RCPP
}
// R_matrix_normalize_sparse
Rcpp::List R_matrix_normalize_sparse(arma::sp_mat mat_norm);
RcppExport SEXP _inferCSN_R_matrix_normalize_sparse(SEXP mat_normSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::sp_mat >::type mat_norm(mat_normSEXP);
    rcpp_result_gen = Rcpp::wrap(R_matrix_normalize_sparse(mat_norm));
    return rcpp_result_gen;
END_RCPP
}
// R_matrix_center_dense
Rcpp::List R_matrix_center_dense(const arma::mat mat, arma::mat X_normalized, bool intercept);
RcppExport SEXP _inferCSN_R_matrix_center_dense(SEXP matSEXP, SEXP X_normalizedSEXP, SEXP interceptSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat >::type mat(matSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X_normalized(X_normalizedSEXP);
    Rcpp::traits::input_parameter< bool >::type intercept(interceptSEXP);
    rcpp_result_gen = Rcpp::wrap(R_matrix_center_dense(mat, X_normalized, intercept));
    return rcpp_result_gen;
END_RCPP
}
// R_matrix_center_sparse
Rcpp::List R_matrix_center_sparse(const arma::sp_mat mat, arma::sp_mat X_normalized, bool intercept);
RcppExport SEXP _inferCSN_R_matrix_center_sparse(SEXP matSEXP, SEXP X_normalizedSEXP, SEXP interceptSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::sp_mat >::type mat(matSEXP);
    Rcpp::traits::input_parameter< arma::sp_mat >::type X_normalized(X_normalizedSEXP);
    Rcpp::traits::input_parameter< bool >::type intercept(interceptSEXP);
    rcpp_result_gen = Rcpp::wrap(R_matrix_center_sparse(mat, X_normalized, intercept));
    return rcpp_result_gen;
END_RCPP
}
// filter_sort_matrix
NumericMatrix filter_sort_matrix(NumericMatrix network_matrix, Nullable<CharacterVector> regulators, Nullable<CharacterVector> targets);
RcppExport SEXP _inferCSN_filter_sort_matrix(SEXP network_matrixSEXP, SEXP regulatorsSEXP, SEXP targetsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type network_matrix(network_matrixSEXP);
    Rcpp::traits::input_parameter< Nullable<CharacterVector> >::type regulators(regulatorsSEXP);
    Rcpp::traits::input_parameter< Nullable<CharacterVector> >::type targets(targetsSEXP);
    rcpp_result_gen = Rcpp::wrap(filter_sort_matrix(network_matrix, regulators, targets));
    return rcpp_result_gen;
END_RCPP
}
// matrix_to_table
DataFrame matrix_to_table(NumericMatrix network_matrix, Nullable<CharacterVector> regulators, Nullable<CharacterVector> targets, double threshold);
RcppExport SEXP _inferCSN_matrix_to_table(SEXP network_matrixSEXP, SEXP regulatorsSEXP, SEXP targetsSEXP, SEXP thresholdSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type network_matrix(network_matrixSEXP);
    Rcpp::traits::input_parameter< Nullable<CharacterVector> >::type regulators(regulatorsSEXP);
    Rcpp::traits::input_parameter< Nullable<CharacterVector> >::type targets(targetsSEXP);
    Rcpp::traits::input_parameter< double >::type threshold(thresholdSEXP);
    rcpp_result_gen = Rcpp::wrap(matrix_to_table(network_matrix, regulators, targets, threshold));
    return rcpp_result_gen;
END_RCPP
}
// prepare_calculate_metrics
DataFrame prepare_calculate_metrics(DataFrame network_table, DataFrame ground_truth);
RcppExport SEXP _inferCSN_prepare_calculate_metrics(SEXP network_tableSEXP, SEXP ground_truthSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame >::type network_table(network_tableSEXP);
    Rcpp::traits::input_parameter< DataFrame >::type ground_truth(ground_truthSEXP);
    rcpp_result_gen = Rcpp::wrap(prepare_calculate_metrics(network_table, ground_truth));
    return rcpp_result_gen;
END_RCPP
}
// network_format
DataFrame network_format(DataFrame network_table, Nullable<CharacterVector> regulators, Nullable<CharacterVector> targets, bool abs_weight);
RcppExport SEXP _inferCSN_network_format(SEXP network_tableSEXP, SEXP regulatorsSEXP, SEXP targetsSEXP, SEXP abs_weightSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame >::type network_table(network_tableSEXP);
    Rcpp::traits::input_parameter< Nullable<CharacterVector> >::type regulators(regulatorsSEXP);
    Rcpp::traits::input_parameter< Nullable<CharacterVector> >::type targets(targetsSEXP);
    Rcpp::traits::input_parameter< bool >::type abs_weight(abs_weightSEXP);
    rcpp_result_gen = Rcpp::wrap(network_format(network_table, regulators, targets, abs_weight));
    return rcpp_result_gen;
END_RCPP
}
// table_to_matrix
NumericMatrix table_to_matrix(DataFrame network_table, Nullable<CharacterVector> regulators, Nullable<CharacterVector> targets);
RcppExport SEXP _inferCSN_table_to_matrix(SEXP network_tableSEXP, SEXP regulatorsSEXP, SEXP targetsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame >::type network_table(network_tableSEXP);
    Rcpp::traits::input_parameter< Nullable<CharacterVector> >::type regulators(regulatorsSEXP);
    Rcpp::traits::input_parameter< Nullable<CharacterVector> >::type targets(targetsSEXP);
    rcpp_result_gen = Rcpp::wrap(table_to_matrix(network_table, regulators, targets));
    return rcpp_result_gen;
END_RCPP
}
// weight_sift
Rcpp::DataFrame weight_sift(Rcpp::DataFrame table);
RcppExport SEXP _inferCSN_weight_sift(SEXP tableSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::DataFrame >::type table(tableSEXP);
    rcpp_result_gen = Rcpp::wrap(weight_sift(table));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_inferCSN_srm_model_sparse", (DL_FUNC) &_inferCSN_srm_model_sparse, 26},
    {"_inferCSN_srm_model_dense", (DL_FUNC) &_inferCSN_srm_model_dense, 26},
    {"_inferCSN_srm_model_cv_sparse", (DL_FUNC) &_inferCSN_srm_model_cv_sparse, 28},
    {"_inferCSN_srm_model_cv_dense", (DL_FUNC) &_inferCSN_srm_model_cv_dense, 28},
    {"_inferCSN_R_matrix_column_get_dense", (DL_FUNC) &_inferCSN_R_matrix_column_get_dense, 2},
    {"_inferCSN_R_matrix_column_get_sparse", (DL_FUNC) &_inferCSN_R_matrix_column_get_sparse, 2},
    {"_inferCSN_R_matrix_rows_get_dense", (DL_FUNC) &_inferCSN_R_matrix_rows_get_dense, 2},
    {"_inferCSN_R_matrix_rows_get_sparse", (DL_FUNC) &_inferCSN_R_matrix_rows_get_sparse, 2},
    {"_inferCSN_R_matrix_vector_schur_product_dense", (DL_FUNC) &_inferCSN_R_matrix_vector_schur_product_dense, 2},
    {"_inferCSN_R_matrix_vector_schur_product_sparse", (DL_FUNC) &_inferCSN_R_matrix_vector_schur_product_sparse, 2},
    {"_inferCSN_R_matrix_vector_divide_dense", (DL_FUNC) &_inferCSN_R_matrix_vector_divide_dense, 2},
    {"_inferCSN_R_matrix_vector_divide_sparse", (DL_FUNC) &_inferCSN_R_matrix_vector_divide_sparse, 2},
    {"_inferCSN_R_matrix_column_sums_dense", (DL_FUNC) &_inferCSN_R_matrix_column_sums_dense, 1},
    {"_inferCSN_R_matrix_column_sums_sparse", (DL_FUNC) &_inferCSN_R_matrix_column_sums_sparse, 1},
    {"_inferCSN_R_matrix_column_dot_dense", (DL_FUNC) &_inferCSN_R_matrix_column_dot_dense, 3},
    {"_inferCSN_R_matrix_column_dot_sparse", (DL_FUNC) &_inferCSN_R_matrix_column_dot_sparse, 3},
    {"_inferCSN_R_matrix_column_mult_dense", (DL_FUNC) &_inferCSN_R_matrix_column_mult_dense, 3},
    {"_inferCSN_R_matrix_column_mult_sparse", (DL_FUNC) &_inferCSN_R_matrix_column_mult_sparse, 3},
    {"_inferCSN_R_matrix_normalize_dense", (DL_FUNC) &_inferCSN_R_matrix_normalize_dense, 1},
    {"_inferCSN_R_matrix_normalize_sparse", (DL_FUNC) &_inferCSN_R_matrix_normalize_sparse, 1},
    {"_inferCSN_R_matrix_center_dense", (DL_FUNC) &_inferCSN_R_matrix_center_dense, 3},
    {"_inferCSN_R_matrix_center_sparse", (DL_FUNC) &_inferCSN_R_matrix_center_sparse, 3},
    {"_inferCSN_filter_sort_matrix", (DL_FUNC) &_inferCSN_filter_sort_matrix, 3},
    {"_inferCSN_matrix_to_table", (DL_FUNC) &_inferCSN_matrix_to_table, 4},
    {"_inferCSN_prepare_calculate_metrics", (DL_FUNC) &_inferCSN_prepare_calculate_metrics, 2},
    {"_inferCSN_network_format", (DL_FUNC) &_inferCSN_network_format, 4},
    {"_inferCSN_table_to_matrix", (DL_FUNC) &_inferCSN_table_to_matrix, 3},
    {"_inferCSN_weight_sift", (DL_FUNC) &_inferCSN_weight_sift, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_inferCSN(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
