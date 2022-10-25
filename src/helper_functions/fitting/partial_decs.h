// Additional declarations for the partial derivative functions

#ifndef PARTIAL_DECS_H
#define PARTIAL_DECS_H

#include <Rcpp.h>
// [[Rcpp::plugins(cpp11)]]
#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

#include "../parameter_checks.h"

using std::isnormal;
using Rcpp::NumericVector;

//--------- Templates for Partial Derivatives --------------------------------//
typedef double (*ParFunc)(const double&, const double&, const double&,
                          const double&, const double&, const double&,
                          const double&);

//--------------------------- Function Declarations --------------------------//
vector<double> partial_pdf(const ParFunc& parf,
                           const NumericVector& rt,
                           const SEXP& response,
                           const NumericVector& v,
                           const NumericVector& a,
                           const NumericVector& t0,
                           const NumericVector& w,
                           const NumericVector& sv,
                           const NumericVector& sigma,
                           const double& sl_thresh,
                           NumericVector err_tol);

#endif // PARTIAL_DECS_H
