// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// cpp_dfddm_fast
NumericVector cpp_dfddm_fast(const NumericVector& rt, const double& a, const double& v, const double& t0, const double& w, const double& eps);
RcppExport SEXP _fddm_cpp_dfddm_fast(SEXP rtSEXP, SEXP aSEXP, SEXP vSEXP, SEXP t0SEXP, SEXP wSEXP, SEXP epsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type rt(rtSEXP);
    Rcpp::traits::input_parameter< const double& >::type a(aSEXP);
    Rcpp::traits::input_parameter< const double& >::type v(vSEXP);
    Rcpp::traits::input_parameter< const double& >::type t0(t0SEXP);
    Rcpp::traits::input_parameter< const double& >::type w(wSEXP);
    Rcpp::traits::input_parameter< const double& >::type eps(epsSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_dfddm_fast(rt, a, v, t0, w, eps));
    return rcpp_result_gen;
END_RCPP
}
// cpp_dfddm
NumericVector cpp_dfddm(const NumericVector& rt, const LogicalVector& response, const NumericVector& a, const NumericVector& v, const NumericVector& t0, const NumericVector& w, const NumericVector& sv, const LogicalVector& log_prob, const std::string& n_terms_small, const std::string& summation_small, const std::string& scale, const NumericVector& eps);
RcppExport SEXP _fddm_cpp_dfddm(SEXP rtSEXP, SEXP responseSEXP, SEXP aSEXP, SEXP vSEXP, SEXP t0SEXP, SEXP wSEXP, SEXP svSEXP, SEXP log_probSEXP, SEXP n_terms_smallSEXP, SEXP summation_smallSEXP, SEXP scaleSEXP, SEXP epsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type rt(rtSEXP);
    Rcpp::traits::input_parameter< const LogicalVector& >::type response(responseSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type a(aSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type v(vSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type t0(t0SEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type w(wSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type sv(svSEXP);
    Rcpp::traits::input_parameter< const LogicalVector& >::type log_prob(log_probSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type n_terms_small(n_terms_smallSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type summation_small(summation_smallSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type scale(scaleSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type eps(epsSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_dfddm(rt, response, a, v, t0, w, sv, log_prob, n_terms_small, summation_small, scale, eps));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_fddm_cpp_dfddm_fast", (DL_FUNC) &_fddm_cpp_dfddm_fast, 6},
    {"_fddm_cpp_dfddm", (DL_FUNC) &_fddm_cpp_dfddm, 12},
    {NULL, NULL, 0}
};

RcppExport void R_init_fddm(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
