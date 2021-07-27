// Header file to declare constants and functions
// This file must have a .h extension rather than .hpp extension because then...
// ...R CMD Check will fail

#include <Rcpp.h>
// [[Rcpp::plugins(cpp11)]]
#include <cmath>
#include <Rmath.h>

using std::vector;
using std::log;
using std::exp;
using std::erf;
using std::sqrt;
using std::isfinite;
using std::isnan;
using std::isnormal;
using Rcpp::NumericVector;



// Constants

// static const double SV_THRESH = 0; // threshold for using variable drift rate
static const double LOG_1_2 = log(0.5);
static const double LOG_100 = log(100);
static const double PI_CONST = 3.14159265358979323846; // define pi like C++
static const double SQRT_2PI = sqrt(2 * PI_CONST);
static const double SQRT_2PI_INV = 1 / SQRT_2PI;
static const double SQRT_2_INV_NEG = -1 / sqrt(2);
static const char EMPTYCHAR = '\0'; // literally, the empty character



// define types for functions in if-else in pfddm.cpp
typedef double (*ErfFunc)(const double& arg);
typedef double (*MillsFunc)(const double& x);
typedef double (*DisFunc)(const double& t, const double& a, const double& v,
                          const double& w, const double& sv, const double& err);



// Mills Ratio (and Approximation)
double r_mills(const double& x);
double c_mills(const double& x);
double zeta_mills(const double& x);



// Infinite Summation Approximations
double mills_sum(const double& t, const double& a, const double& v,
                 const double& w, const double& sv, const double& err);
double ncdf_sum(const double& t, const double& a, const double& v,
                const double& w, const double& sv, const double& err);



// Density Functions
double Fs_mills(const double& t, const double& a, const double& v,
                const double& w, const double& sv, const double& err);
double Fs_mills_log(const double& t, const double& a, const double& v,
                    const double& w, const double& sv, const double& err);
double Fs_ncdf(const double& t, const double& a, const double& v,
               const double& w, const double& sv, const double& err);
double Fs_ncdf_log(const double& t, const double& a, const double& v,
                   const double& w, const double& sv, const double& err);



// Helper Functions
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
void determine_method(const string& method, DisFunc& disf,
                      double& rt0, const bool& log_prob);
double prob_lower(const double& a, const double& v, const double& w,
                  const double& rt0);
void calculate_cdf(const int& Nrt, const int& Na, const int& Nv, const int& Nt0,
                   const int& Nw, const int& Nsv, const int& Nsig,
                   const int& Nerr, const int& Nmax,
                   const NumericVector& rt,
                   const NumericVector& a, const NumericVector& v,
                   const NumericVector& t0, const NumericVector& w,
                   const NumericVector& sv, const NumericVector& sigma,
                   const NumericVector& err, vector<double>& out,
                   const double& rt0, const DisFunc& disf);
