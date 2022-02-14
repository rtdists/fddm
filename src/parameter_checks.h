// Functions to check valid parameter values in dfddm and pfddm

#ifndef PARAMETER_CHECKS_H
#define PARAMETER_CHECKS_H

#include <Rcpp.h>
// [[Rcpp::plugins(cpp11)]]
#include <cmath>

using std::vector;
using std::string;
using std::max;
using std::to_string;
using std::isfinite;
using std::isnan;
using Rcpp::stop;
using Rcpp::warning;
using Rcpp::NumericVector;
using Rcpp::LogicalVector;


void convert_responses(const SEXP& response, int& Nres, int& Nmax,
                       vector<double>& out, const double& rt0, bool& valid);

bool parameter_check(const int& Nrt, int& Nres, const int& Na, const int& Nv,
                     const int& Nt0, const int& Nw, const int& Nsv,
                     const int& Nsig, const int& Nerr, int& Nmax,
                     const NumericVector& rt, const SEXP& response,
                     const NumericVector& a, const NumericVector& v,
                     const NumericVector& t0, const NumericVector& w,
                     const NumericVector& sv, const NumericVector& sigma,
                     const NumericVector& err,
                     vector<double>& out, const double& rt0);

#endif // PARAMETER_CHECKS_H
