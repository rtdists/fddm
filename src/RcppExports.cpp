// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// dt_dfddm
NumericVector dt_dfddm(const NumericVector& rt, const SEXP& response, const NumericVector& a, const NumericVector& v, const NumericVector& t0, const NumericVector& w, const NumericVector& sv, const NumericVector& sigma, const double& sl_thresh, const NumericVector& err_tol);
RcppExport SEXP _fddm_dt_dfddm(SEXP rtSEXP, SEXP responseSEXP, SEXP aSEXP, SEXP vSEXP, SEXP t0SEXP, SEXP wSEXP, SEXP svSEXP, SEXP sigmaSEXP, SEXP sl_threshSEXP, SEXP err_tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type rt(rtSEXP);
    Rcpp::traits::input_parameter< const SEXP& >::type response(responseSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type a(aSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type v(vSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type t0(t0SEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type w(wSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type sv(svSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< const double& >::type sl_thresh(sl_threshSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type err_tol(err_tolSEXP);
    rcpp_result_gen = Rcpp::wrap(dt_dfddm(rt, response, a, v, t0, w, sv, sigma, sl_thresh, err_tol));
    return rcpp_result_gen;
END_RCPP
}
// dt0_dfddm
NumericVector dt0_dfddm(const NumericVector& rt, const SEXP& response, const NumericVector& a, const NumericVector& v, const NumericVector& t0, const NumericVector& w, const NumericVector& sv, const NumericVector& sigma, const double& sl_thresh, const NumericVector& err_tol);
RcppExport SEXP _fddm_dt0_dfddm(SEXP rtSEXP, SEXP responseSEXP, SEXP aSEXP, SEXP vSEXP, SEXP t0SEXP, SEXP wSEXP, SEXP svSEXP, SEXP sigmaSEXP, SEXP sl_threshSEXP, SEXP err_tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type rt(rtSEXP);
    Rcpp::traits::input_parameter< const SEXP& >::type response(responseSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type a(aSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type v(vSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type t0(t0SEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type w(wSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type sv(svSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< const double& >::type sl_thresh(sl_threshSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type err_tol(err_tolSEXP);
    rcpp_result_gen = Rcpp::wrap(dt0_dfddm(rt, response, a, v, t0, w, sv, sigma, sl_thresh, err_tol));
    return rcpp_result_gen;
END_RCPP
}
// da_dfddm
NumericVector da_dfddm(const NumericVector& rt, const SEXP& response, const NumericVector& a, const NumericVector& v, const NumericVector& t0, const NumericVector& w, const NumericVector& sv, const NumericVector& sigma, const double& sl_thresh, const NumericVector& err_tol);
RcppExport SEXP _fddm_da_dfddm(SEXP rtSEXP, SEXP responseSEXP, SEXP aSEXP, SEXP vSEXP, SEXP t0SEXP, SEXP wSEXP, SEXP svSEXP, SEXP sigmaSEXP, SEXP sl_threshSEXP, SEXP err_tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type rt(rtSEXP);
    Rcpp::traits::input_parameter< const SEXP& >::type response(responseSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type a(aSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type v(vSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type t0(t0SEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type w(wSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type sv(svSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< const double& >::type sl_thresh(sl_threshSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type err_tol(err_tolSEXP);
    rcpp_result_gen = Rcpp::wrap(da_dfddm(rt, response, a, v, t0, w, sv, sigma, sl_thresh, err_tol));
    return rcpp_result_gen;
END_RCPP
}
// dv_dfddm
NumericVector dv_dfddm(const NumericVector& rt, const SEXP& response, const NumericVector& a, const NumericVector& v, const NumericVector& t0, const NumericVector& w, const NumericVector& sv, const NumericVector& sigma, const double& sl_thresh, const NumericVector& err_tol);
RcppExport SEXP _fddm_dv_dfddm(SEXP rtSEXP, SEXP responseSEXP, SEXP aSEXP, SEXP vSEXP, SEXP t0SEXP, SEXP wSEXP, SEXP svSEXP, SEXP sigmaSEXP, SEXP sl_threshSEXP, SEXP err_tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type rt(rtSEXP);
    Rcpp::traits::input_parameter< const SEXP& >::type response(responseSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type a(aSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type v(vSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type t0(t0SEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type w(wSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type sv(svSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< const double& >::type sl_thresh(sl_threshSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type err_tol(err_tolSEXP);
    rcpp_result_gen = Rcpp::wrap(dv_dfddm(rt, response, a, v, t0, w, sv, sigma, sl_thresh, err_tol));
    return rcpp_result_gen;
END_RCPP
}
// dw_dfddm
NumericVector dw_dfddm(const NumericVector& rt, const SEXP& response, const NumericVector& a, const NumericVector& v, const NumericVector& t0, const NumericVector& w, const NumericVector& sv, const NumericVector& sigma, const double& sl_thresh, const NumericVector& err_tol);
RcppExport SEXP _fddm_dw_dfddm(SEXP rtSEXP, SEXP responseSEXP, SEXP aSEXP, SEXP vSEXP, SEXP t0SEXP, SEXP wSEXP, SEXP svSEXP, SEXP sigmaSEXP, SEXP sl_threshSEXP, SEXP err_tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type rt(rtSEXP);
    Rcpp::traits::input_parameter< const SEXP& >::type response(responseSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type a(aSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type v(vSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type t0(t0SEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type w(wSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type sv(svSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< const double& >::type sl_thresh(sl_threshSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type err_tol(err_tolSEXP);
    rcpp_result_gen = Rcpp::wrap(dw_dfddm(rt, response, a, v, t0, w, sv, sigma, sl_thresh, err_tol));
    return rcpp_result_gen;
END_RCPP
}
// dsv_dfddm
NumericVector dsv_dfddm(const NumericVector& rt, const SEXP& response, const NumericVector& a, const NumericVector& v, const NumericVector& t0, const NumericVector& w, const NumericVector& sv, const NumericVector& sigma, const double& sl_thresh, const NumericVector& err_tol);
RcppExport SEXP _fddm_dsv_dfddm(SEXP rtSEXP, SEXP responseSEXP, SEXP aSEXP, SEXP vSEXP, SEXP t0SEXP, SEXP wSEXP, SEXP svSEXP, SEXP sigmaSEXP, SEXP sl_threshSEXP, SEXP err_tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type rt(rtSEXP);
    Rcpp::traits::input_parameter< const SEXP& >::type response(responseSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type a(aSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type v(vSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type t0(t0SEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type w(wSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type sv(svSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< const double& >::type sl_thresh(sl_threshSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type err_tol(err_tolSEXP);
    rcpp_result_gen = Rcpp::wrap(dsv_dfddm(rt, response, a, v, t0, w, sv, sigma, sl_thresh, err_tol));
    return rcpp_result_gen;
END_RCPP
}
// dfddm
NumericVector dfddm(const NumericVector& rt, const SEXP& response, const NumericVector& a, const NumericVector& v, const NumericVector& t0, const NumericVector& w, const NumericVector& sv, const NumericVector& sigma, const NumericVector& err_tol, const bool& log, const std::string& switch_mech, double switch_thresh, const std::string& n_terms_small, const std::string& summation_small);
RcppExport SEXP _fddm_dfddm(SEXP rtSEXP, SEXP responseSEXP, SEXP aSEXP, SEXP vSEXP, SEXP t0SEXP, SEXP wSEXP, SEXP svSEXP, SEXP sigmaSEXP, SEXP err_tolSEXP, SEXP logSEXP, SEXP switch_mechSEXP, SEXP switch_threshSEXP, SEXP n_terms_smallSEXP, SEXP summation_smallSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type rt(rtSEXP);
    Rcpp::traits::input_parameter< const SEXP& >::type response(responseSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type a(aSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type v(vSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type t0(t0SEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type w(wSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type sv(svSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type err_tol(err_tolSEXP);
    Rcpp::traits::input_parameter< const bool& >::type log(logSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type switch_mech(switch_mechSEXP);
    Rcpp::traits::input_parameter< double >::type switch_thresh(switch_threshSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type n_terms_small(n_terms_smallSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type summation_small(summation_smallSEXP);
    rcpp_result_gen = Rcpp::wrap(dfddm(rt, response, a, v, t0, w, sv, sigma, err_tol, log, switch_mech, switch_thresh, n_terms_small, summation_small));
    return rcpp_result_gen;
END_RCPP
}
// pfddm
NumericVector pfddm(const NumericVector& rt, const SEXP& response, const NumericVector& a, const NumericVector& v, const NumericVector& t0, const NumericVector& w, const NumericVector& sv, const NumericVector& sigma, NumericVector err_tol, const bool& log, const std::string& method);
RcppExport SEXP _fddm_pfddm(SEXP rtSEXP, SEXP responseSEXP, SEXP aSEXP, SEXP vSEXP, SEXP t0SEXP, SEXP wSEXP, SEXP svSEXP, SEXP sigmaSEXP, SEXP err_tolSEXP, SEXP logSEXP, SEXP methodSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type rt(rtSEXP);
    Rcpp::traits::input_parameter< const SEXP& >::type response(responseSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type a(aSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type v(vSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type t0(t0SEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type w(wSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type sv(svSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type err_tol(err_tolSEXP);
    Rcpp::traits::input_parameter< const bool& >::type log(logSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type method(methodSEXP);
    rcpp_result_gen = Rcpp::wrap(pfddm(rt, response, a, v, t0, w, sv, sigma, err_tol, log, method));
    return rcpp_result_gen;
END_RCPP
}

RcppExport SEXP _rcpp_module_boot_fddm_fit();

static const R_CallMethodDef CallEntries[] = {
    {"_fddm_dt_dfddm", (DL_FUNC) &_fddm_dt_dfddm, 10},
    {"_fddm_dt0_dfddm", (DL_FUNC) &_fddm_dt0_dfddm, 10},
    {"_fddm_da_dfddm", (DL_FUNC) &_fddm_da_dfddm, 10},
    {"_fddm_dv_dfddm", (DL_FUNC) &_fddm_dv_dfddm, 10},
    {"_fddm_dw_dfddm", (DL_FUNC) &_fddm_dw_dfddm, 10},
    {"_fddm_dsv_dfddm", (DL_FUNC) &_fddm_dsv_dfddm, 10},
    {"_fddm_dfddm", (DL_FUNC) &_fddm_dfddm, 14},
    {"_fddm_pfddm", (DL_FUNC) &_fddm_pfddm, 11},
    {"_rcpp_module_boot_fddm_fit", (DL_FUNC) &_rcpp_module_boot_fddm_fit, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_fddm(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
