// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// RunModularityClusteringCpp
IntegerVector RunModularityClusteringCpp(Eigen::SparseMatrix<double> SNN, int modularityFunction, double resolution, int algorithm, int nRandomStarts, int nIterations, int randomSeed, bool printOutput);
RcppExport SEXP _scMRMA_RunModularityClusteringCpp(SEXP SNNSEXP, SEXP modularityFunctionSEXP, SEXP resolutionSEXP, SEXP algorithmSEXP, SEXP nRandomStartsSEXP, SEXP nIterationsSEXP, SEXP randomSeedSEXP, SEXP printOutputSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::SparseMatrix<double> >::type SNN(SNNSEXP);
    Rcpp::traits::input_parameter< int >::type modularityFunction(modularityFunctionSEXP);
    Rcpp::traits::input_parameter< double >::type resolution(resolutionSEXP);
    Rcpp::traits::input_parameter< int >::type algorithm(algorithmSEXP);
    Rcpp::traits::input_parameter< int >::type nRandomStarts(nRandomStartsSEXP);
    Rcpp::traits::input_parameter< int >::type nIterations(nIterationsSEXP);
    Rcpp::traits::input_parameter< int >::type randomSeed(randomSeedSEXP);
    Rcpp::traits::input_parameter< bool >::type printOutput(printOutputSEXP);
    rcpp_result_gen = Rcpp::wrap(RunModularityClusteringCpp(SNN, modularityFunction, resolution, algorithm, nRandomStarts, nIterations, randomSeed, printOutput));
    return rcpp_result_gen;
END_RCPP
}
// getSNN
Eigen::SparseMatrix<double> getSNN(Eigen::MatrixXd knnMatrix, double prune);
RcppExport SEXP _scMRMA_getSNN(SEXP knnMatrixSEXP, SEXP pruneSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type knnMatrix(knnMatrixSEXP);
    Rcpp::traits::input_parameter< double >::type prune(pruneSEXP);
    rcpp_result_gen = Rcpp::wrap(getSNN(knnMatrix, prune));
    return rcpp_result_gen;
END_RCPP
}
// LogNorm
Eigen::SparseMatrix<double> LogNorm(Eigen::SparseMatrix<double> data, bool display_progress);
RcppExport SEXP _scMRMA_LogNorm(SEXP dataSEXP, SEXP display_progressSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::SparseMatrix<double> >::type data(dataSEXP);
    Rcpp::traits::input_parameter< bool >::type display_progress(display_progressSEXP);
    rcpp_result_gen = Rcpp::wrap(LogNorm(data, display_progress));
    return rcpp_result_gen;
END_RCPP
}
// SparseRowVar2
NumericVector SparseRowVar2(Eigen::SparseMatrix<double> mat, NumericVector mu, bool display_progress);
RcppExport SEXP _scMRMA_SparseRowVar2(SEXP matSEXP, SEXP muSEXP, SEXP display_progressSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::SparseMatrix<double> >::type mat(matSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< bool >::type display_progress(display_progressSEXP);
    rcpp_result_gen = Rcpp::wrap(SparseRowVar2(mat, mu, display_progress));
    return rcpp_result_gen;
END_RCPP
}
// SparseRowVarStd
NumericVector SparseRowVarStd(Eigen::SparseMatrix<double> mat, NumericVector mu, NumericVector sd, double vmax, bool display_progress);
RcppExport SEXP _scMRMA_SparseRowVarStd(SEXP matSEXP, SEXP muSEXP, SEXP sdSEXP, SEXP vmaxSEXP, SEXP display_progressSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::SparseMatrix<double> >::type mat(matSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sd(sdSEXP);
    Rcpp::traits::input_parameter< double >::type vmax(vmaxSEXP);
    Rcpp::traits::input_parameter< bool >::type display_progress(display_progressSEXP);
    rcpp_result_gen = Rcpp::wrap(SparseRowVarStd(mat, mu, sd, vmax, display_progress));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_scMRMA_RunModularityClusteringCpp", (DL_FUNC) &_scMRMA_RunModularityClusteringCpp, 8},
    {"_scMRMA_getSNN", (DL_FUNC) &_scMRMA_getSNN, 2},
    {"_scMRMA_LogNorm", (DL_FUNC) &_scMRMA_LogNorm, 2},
    {"_scMRMA_SparseRowVar2", (DL_FUNC) &_scMRMA_SparseRowVar2, 3},
    {"_scMRMA_SparseRowVarStd", (DL_FUNC) &_scMRMA_SparseRowVarStd, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_scMRMA(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
