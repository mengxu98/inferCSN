// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// SRM_model_fit_sparse
Rcpp::List SRM_model_fit_sparse(const arma::sp_mat& X, const arma::vec& y, const std::string Loss, const std::string Penalty, const std::string Algorithm, const std::size_t NnzStopNum, const std::size_t G_ncols, const std::size_t G_nrows, const double Lambda2Max, const double Lambda2Min, const bool PartialSort, const std::size_t MaxIters, const double rtol, const double atol, const bool ActiveSet, const std::size_t ActiveSetNum, const std::size_t MaxNumSwaps, const double ScaleDownFactor, const std::size_t ScreenSize, const bool LambdaU, const std::vector< std::vector<double> > Lambdas, const std::size_t ExcludeFirstK, const bool Intercept, const bool withBounds, const arma::vec& Lows, const arma::vec& Highs);
RcppExport SEXP _inferCSN_SRM_model_fit_sparse(SEXP XSEXP, SEXP ySEXP, SEXP LossSEXP, SEXP PenaltySEXP, SEXP AlgorithmSEXP, SEXP NnzStopNumSEXP, SEXP G_ncolsSEXP, SEXP G_nrowsSEXP, SEXP Lambda2MaxSEXP, SEXP Lambda2MinSEXP, SEXP PartialSortSEXP, SEXP MaxItersSEXP, SEXP rtolSEXP, SEXP atolSEXP, SEXP ActiveSetSEXP, SEXP ActiveSetNumSEXP, SEXP MaxNumSwapsSEXP, SEXP ScaleDownFactorSEXP, SEXP ScreenSizeSEXP, SEXP LambdaUSEXP, SEXP LambdasSEXP, SEXP ExcludeFirstKSEXP, SEXP InterceptSEXP, SEXP withBoundsSEXP, SEXP LowsSEXP, SEXP HighsSEXP) {
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
    Rcpp::traits::input_parameter< const std::vector< std::vector<double> > >::type Lambdas(LambdasSEXP);
    Rcpp::traits::input_parameter< const std::size_t >::type ExcludeFirstK(ExcludeFirstKSEXP);
    Rcpp::traits::input_parameter< const bool >::type Intercept(InterceptSEXP);
    Rcpp::traits::input_parameter< const bool >::type withBounds(withBoundsSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type Lows(LowsSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type Highs(HighsSEXP);
    rcpp_result_gen = Rcpp::wrap(SRM_model_fit_sparse(X, y, Loss, Penalty, Algorithm, NnzStopNum, G_ncols, G_nrows, Lambda2Max, Lambda2Min, PartialSort, MaxIters, rtol, atol, ActiveSet, ActiveSetNum, MaxNumSwaps, ScaleDownFactor, ScreenSize, LambdaU, Lambdas, ExcludeFirstK, Intercept, withBounds, Lows, Highs));
    return rcpp_result_gen;
END_RCPP
}
// SRM_model_fit_dense
Rcpp::List SRM_model_fit_dense(const arma::mat& X, const arma::vec& y, const std::string Loss, const std::string Penalty, const std::string Algorithm, const std::size_t NnzStopNum, const std::size_t G_ncols, const std::size_t G_nrows, const double Lambda2Max, const double Lambda2Min, const bool PartialSort, const std::size_t MaxIters, const double rtol, const double atol, const bool ActiveSet, const std::size_t ActiveSetNum, const std::size_t MaxNumSwaps, const double ScaleDownFactor, const std::size_t ScreenSize, const bool LambdaU, const std::vector< std::vector<double> > Lambdas, const std::size_t ExcludeFirstK, const bool Intercept, const bool withBounds, const arma::vec& Lows, const arma::vec& Highs);
RcppExport SEXP _inferCSN_SRM_model_fit_dense(SEXP XSEXP, SEXP ySEXP, SEXP LossSEXP, SEXP PenaltySEXP, SEXP AlgorithmSEXP, SEXP NnzStopNumSEXP, SEXP G_ncolsSEXP, SEXP G_nrowsSEXP, SEXP Lambda2MaxSEXP, SEXP Lambda2MinSEXP, SEXP PartialSortSEXP, SEXP MaxItersSEXP, SEXP rtolSEXP, SEXP atolSEXP, SEXP ActiveSetSEXP, SEXP ActiveSetNumSEXP, SEXP MaxNumSwapsSEXP, SEXP ScaleDownFactorSEXP, SEXP ScreenSizeSEXP, SEXP LambdaUSEXP, SEXP LambdasSEXP, SEXP ExcludeFirstKSEXP, SEXP InterceptSEXP, SEXP withBoundsSEXP, SEXP LowsSEXP, SEXP HighsSEXP) {
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
    Rcpp::traits::input_parameter< const std::vector< std::vector<double> > >::type Lambdas(LambdasSEXP);
    Rcpp::traits::input_parameter< const std::size_t >::type ExcludeFirstK(ExcludeFirstKSEXP);
    Rcpp::traits::input_parameter< const bool >::type Intercept(InterceptSEXP);
    Rcpp::traits::input_parameter< const bool >::type withBounds(withBoundsSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type Lows(LowsSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type Highs(HighsSEXP);
    rcpp_result_gen = Rcpp::wrap(SRM_model_fit_dense(X, y, Loss, Penalty, Algorithm, NnzStopNum, G_ncols, G_nrows, Lambda2Max, Lambda2Min, PartialSort, MaxIters, rtol, atol, ActiveSet, ActiveSetNum, MaxNumSwaps, ScaleDownFactor, ScreenSize, LambdaU, Lambdas, ExcludeFirstK, Intercept, withBounds, Lows, Highs));
    return rcpp_result_gen;
END_RCPP
}
// SRM_model_fit_CV_sparse
Rcpp::List SRM_model_fit_CV_sparse(const arma::sp_mat& X, const arma::vec& y, const std::string Loss, const std::string Penalty, const std::string Algorithm, const std::size_t NnzStopNum, const std::size_t G_ncols, const std::size_t G_nrows, const double Lambda2Max, const double Lambda2Min, const bool PartialSort, const std::size_t MaxIters, const double rtol, const double atol, const bool ActiveSet, const std::size_t ActiveSetNum, const std::size_t MaxNumSwaps, const double ScaleDownFactor, const std::size_t ScreenSize, const bool LambdaU, const std::vector< std::vector<double> > Lambdas, const std::size_t nfolds, const double seed, const std::size_t ExcludeFirstK, const bool Intercept, const bool withBounds, const arma::vec& Lows, const arma::vec& Highs);
RcppExport SEXP _inferCSN_SRM_model_fit_CV_sparse(SEXP XSEXP, SEXP ySEXP, SEXP LossSEXP, SEXP PenaltySEXP, SEXP AlgorithmSEXP, SEXP NnzStopNumSEXP, SEXP G_ncolsSEXP, SEXP G_nrowsSEXP, SEXP Lambda2MaxSEXP, SEXP Lambda2MinSEXP, SEXP PartialSortSEXP, SEXP MaxItersSEXP, SEXP rtolSEXP, SEXP atolSEXP, SEXP ActiveSetSEXP, SEXP ActiveSetNumSEXP, SEXP MaxNumSwapsSEXP, SEXP ScaleDownFactorSEXP, SEXP ScreenSizeSEXP, SEXP LambdaUSEXP, SEXP LambdasSEXP, SEXP nfoldsSEXP, SEXP seedSEXP, SEXP ExcludeFirstKSEXP, SEXP InterceptSEXP, SEXP withBoundsSEXP, SEXP LowsSEXP, SEXP HighsSEXP) {
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
    Rcpp::traits::input_parameter< const std::vector< std::vector<double> > >::type Lambdas(LambdasSEXP);
    Rcpp::traits::input_parameter< const std::size_t >::type nfolds(nfoldsSEXP);
    Rcpp::traits::input_parameter< const double >::type seed(seedSEXP);
    Rcpp::traits::input_parameter< const std::size_t >::type ExcludeFirstK(ExcludeFirstKSEXP);
    Rcpp::traits::input_parameter< const bool >::type Intercept(InterceptSEXP);
    Rcpp::traits::input_parameter< const bool >::type withBounds(withBoundsSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type Lows(LowsSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type Highs(HighsSEXP);
    rcpp_result_gen = Rcpp::wrap(SRM_model_fit_CV_sparse(X, y, Loss, Penalty, Algorithm, NnzStopNum, G_ncols, G_nrows, Lambda2Max, Lambda2Min, PartialSort, MaxIters, rtol, atol, ActiveSet, ActiveSetNum, MaxNumSwaps, ScaleDownFactor, ScreenSize, LambdaU, Lambdas, nfolds, seed, ExcludeFirstK, Intercept, withBounds, Lows, Highs));
    return rcpp_result_gen;
END_RCPP
}
// SRM_model_fit_CV_dense
Rcpp::List SRM_model_fit_CV_dense(const arma::mat& X, const arma::vec& y, const std::string Loss, const std::string Penalty, const std::string Algorithm, const std::size_t NnzStopNum, const std::size_t G_ncols, const std::size_t G_nrows, const double Lambda2Max, const double Lambda2Min, const bool PartialSort, const std::size_t MaxIters, const double rtol, const double atol, const bool ActiveSet, const std::size_t ActiveSetNum, const std::size_t MaxNumSwaps, const double ScaleDownFactor, const std::size_t ScreenSize, const bool LambdaU, const std::vector< std::vector<double> > Lambdas, const std::size_t nfolds, const double seed, const std::size_t ExcludeFirstK, const bool Intercept, const bool withBounds, const arma::vec& Lows, const arma::vec& Highs);
RcppExport SEXP _inferCSN_SRM_model_fit_CV_dense(SEXP XSEXP, SEXP ySEXP, SEXP LossSEXP, SEXP PenaltySEXP, SEXP AlgorithmSEXP, SEXP NnzStopNumSEXP, SEXP G_ncolsSEXP, SEXP G_nrowsSEXP, SEXP Lambda2MaxSEXP, SEXP Lambda2MinSEXP, SEXP PartialSortSEXP, SEXP MaxItersSEXP, SEXP rtolSEXP, SEXP atolSEXP, SEXP ActiveSetSEXP, SEXP ActiveSetNumSEXP, SEXP MaxNumSwapsSEXP, SEXP ScaleDownFactorSEXP, SEXP ScreenSizeSEXP, SEXP LambdaUSEXP, SEXP LambdasSEXP, SEXP nfoldsSEXP, SEXP seedSEXP, SEXP ExcludeFirstKSEXP, SEXP InterceptSEXP, SEXP withBoundsSEXP, SEXP LowsSEXP, SEXP HighsSEXP) {
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
    Rcpp::traits::input_parameter< const std::vector< std::vector<double> > >::type Lambdas(LambdasSEXP);
    Rcpp::traits::input_parameter< const std::size_t >::type nfolds(nfoldsSEXP);
    Rcpp::traits::input_parameter< const double >::type seed(seedSEXP);
    Rcpp::traits::input_parameter< const std::size_t >::type ExcludeFirstK(ExcludeFirstKSEXP);
    Rcpp::traits::input_parameter< const bool >::type Intercept(InterceptSEXP);
    Rcpp::traits::input_parameter< const bool >::type withBounds(withBoundsSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type Lows(LowsSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type Highs(HighsSEXP);
    rcpp_result_gen = Rcpp::wrap(SRM_model_fit_CV_dense(X, y, Loss, Penalty, Algorithm, NnzStopNum, G_ncols, G_nrows, Lambda2Max, Lambda2Min, PartialSort, MaxIters, rtol, atol, ActiveSet, ActiveSetNum, MaxNumSwaps, ScaleDownFactor, ScreenSize, LambdaU, Lambdas, nfolds, seed, ExcludeFirstK, Intercept, withBounds, Lows, Highs));
    return rcpp_result_gen;
END_RCPP
}
// cor_matrix
Rcpp::NumericMatrix cor_matrix(const int p, const double base_cor);
RcppExport SEXP _inferCSN_cor_matrix(SEXP pSEXP, SEXP base_corSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int >::type p(pSEXP);
    Rcpp::traits::input_parameter< const double >::type base_cor(base_corSEXP);
    rcpp_result_gen = Rcpp::wrap(cor_matrix(p, base_cor));
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
// asMatrix
NumericMatrix asMatrix(NumericVector rp, NumericVector cp, NumericVector z, int nrows, int ncols);
RcppExport SEXP _inferCSN_asMatrix(SEXP rpSEXP, SEXP cpSEXP, SEXP zSEXP, SEXP nrowsSEXP, SEXP ncolsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type rp(rpSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type cp(cpSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type z(zSEXP);
    Rcpp::traits::input_parameter< int >::type nrows(nrowsSEXP);
    Rcpp::traits::input_parameter< int >::type ncols(ncolsSEXP);
    rcpp_result_gen = Rcpp::wrap(asMatrix(rp, cp, z, nrows, ncols));
    return rcpp_result_gen;
END_RCPP
}
// asMatrixParallel
NumericMatrix asMatrixParallel(NumericVector rp, NumericVector cp, NumericVector z, int nrows, int ncols);
RcppExport SEXP _inferCSN_asMatrixParallel(SEXP rpSEXP, SEXP cpSEXP, SEXP zSEXP, SEXP nrowsSEXP, SEXP ncolsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type rp(rpSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type cp(cpSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type z(zSEXP);
    Rcpp::traits::input_parameter< int >::type nrows(nrowsSEXP);
    Rcpp::traits::input_parameter< int >::type ncols(ncolsSEXP);
    rcpp_result_gen = Rcpp::wrap(asMatrixParallel(rp, cp, z, nrows, ncols));
    return rcpp_result_gen;
END_RCPP
}
// tableToMatrix
NumericMatrix tableToMatrix(DataFrame weight_table);
RcppExport SEXP _inferCSN_tableToMatrix(SEXP weight_tableSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame >::type weight_table(weight_tableSEXP);
    rcpp_result_gen = Rcpp::wrap(tableToMatrix(weight_table));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_inferCSN_SRM_model_fit_sparse", (DL_FUNC) &_inferCSN_SRM_model_fit_sparse, 26},
    {"_inferCSN_SRM_model_fit_dense", (DL_FUNC) &_inferCSN_SRM_model_fit_dense, 26},
    {"_inferCSN_SRM_model_fit_CV_sparse", (DL_FUNC) &_inferCSN_SRM_model_fit_CV_sparse, 28},
    {"_inferCSN_SRM_model_fit_CV_dense", (DL_FUNC) &_inferCSN_SRM_model_fit_CV_dense, 28},
    {"_inferCSN_cor_matrix", (DL_FUNC) &_inferCSN_cor_matrix, 2},
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
    {"_inferCSN_asMatrix", (DL_FUNC) &_inferCSN_asMatrix, 5},
    {"_inferCSN_asMatrixParallel", (DL_FUNC) &_inferCSN_asMatrixParallel, 5},
    {"_inferCSN_tableToMatrix", (DL_FUNC) &_inferCSN_tableToMatrix, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_inferCSN(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
