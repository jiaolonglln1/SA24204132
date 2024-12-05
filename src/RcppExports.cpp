// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// c_cdmlasso
NumericVector c_cdmlasso(NumericVector y, NumericMatrix z, double lambda, double mu, double tol, int max_iter);
RcppExport SEXP _SA24204132_c_cdmlasso(SEXP ySEXP, SEXP zSEXP, SEXP lambdaSEXP, SEXP muSEXP, SEXP tolSEXP, SEXP max_iterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type z(zSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< int >::type max_iter(max_iterSEXP);
    rcpp_result_gen = Rcpp::wrap(c_cdmlasso(y, z, lambda, mu, tol, max_iter));
    return rcpp_result_gen;
END_RCPP
}
// c_clambda
double c_clambda(NumericVector y, NumericMatrix z, double mu, double learning_rate, int max_iter, double tol);
RcppExport SEXP _SA24204132_c_clambda(SEXP ySEXP, SEXP zSEXP, SEXP muSEXP, SEXP learning_rateSEXP, SEXP max_iterSEXP, SEXP tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type z(zSEXP);
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type learning_rate(learning_rateSEXP);
    Rcpp::traits::input_parameter< int >::type max_iter(max_iterSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    rcpp_result_gen = Rcpp::wrap(c_clambda(y, z, mu, learning_rate, max_iter, tol));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_SA24204132_c_cdmlasso", (DL_FUNC) &_SA24204132_c_cdmlasso, 6},
    {"_SA24204132_c_clambda", (DL_FUNC) &_SA24204132_c_clambda, 6},
    {NULL, NULL, 0}
};

RcppExport void R_init_SA24204132(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}