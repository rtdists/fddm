// Header file to declare constants and functions
// This file must have a .h extension rather than .hpp extension because then
// R CMD Check will fail

#include <Rcpp.h>
// [[Rcpp::plugins(cpp11)]]
#include <cmath>

using std::vector;
using std::string;
using std::log;
using std::exp;
using std::sqrt;
using std::min;
using std::max;
using std::ceil;
using std::isfinite;
using std::isnan;
using std::isnormal;
using Rcpp::stop;
using Rcpp::NumericVector;



// Constants
static const double FOUR_NINTHS = 0.444444444444444444444;
static const double SQRT_2 = 1.41421356237309504880;
static const double SQRT_3 = sqrt(3);
static const double LOG_2 = 0.693147180559945309417;
static const double PI_CONST = 3.14159265358979323846; // define pi like C++
static const double LOG_PI = log(PI_CONST);
static const double O_PI = 0.318309886183790671538;
static const double PI2 = PI_CONST * PI_CONST;
static const double SQRT_2PI = sqrt(2 * PI_CONST);
static const double SQRT_1_2PI = 1 / SQRT_2PI;
static const double SQRT_2_PI = sqrt(2 * O_PI);
static const double SQRT_2_1_PI = SQRT_2 * O_PI;
static const char EMPTYCHAR = '\0'; // literally, the empty character
// INT_MAX is in num_funcs.cpp, maximum value of the int data type = 2147483647



// define types for partial derivative functions
typedef double (*ParFunc)(const double&, const double&, const double&,
                          const double&, const double&, const double&,
                          const double&, const double&);



// Number of Terms
int kl_Nav(const double& taa, const double& err);
int kl_Har(const double& taa, const double& t, const double& err);
int kl_Har_w(const double& taa, const double& t, const double& err);
int ks_Har_w(const double& taa, const double& w, const double& err);



// Infinite Summation Approximations
double sum_small(const double& taa, const double& w, const double& err);
double sum_small_d(const double& taa, const double& w, const double& err);
double sum_small_d_w(const double& taa, const double& w, const int& ks);
double sum_large(const double& taa, const double& w, const int& kl);
double sum_large_d(const double& taa, const double& w, const int& kl);
double sum_large_d_w(const double& taa, const double& w, const int& kl);



// Partial Derivatives (of the PDF) Functions
double pdf_dt(const double& t, const double& resp, const double& a,
              const double& v, const double& w, const double& sv,
              const double& err, const double& sl_thresh);
double pdf_dt0(const double& t, const double& resp, const double& a,
               const double& v, const double& w, const double& sv,
               const double& err, const double& sl_thresh);
double pdf_da(const double& t, const double& resp, const double& a,
              const double& v, const double& w, const double& sv,
              const double& err, const double& sl_thresh);
double pdf_dv(const double& t, const double& resp, const double& a,
              const double& v, const double& w, const double& sv,
              const double& err, const double& sl_thresh);
double pdf_dw(const double& t, const double& resp, const double& a,
              const double& v, const double& w, const double& sv,
              const double& err, const double& sl_thresh);
double pdf_dsv(const double& t, const double& resp, const double& a,
               const double& v, const double& w, const double& sv,
               const double& err, const double& sl_thresh);
// Partial Derivatives of the log(PDF) Functions
// double pdf_dt_log(const double& t, const double& resp, const double& a,
//                   const double& v, const double& w, const double& sv,
//                   const double& err, const double& sl_thresh);
// double pdf_dt0_log(const double& t, const double& resp, const double& a,
//                    const double& v, const double& w, const double& sv,
//                    const double& err, const double& sl_thresh);
// double pdf_da_log(const double& t, const double& resp, const double& a,
//                   const double& v, const double& w, const double& sv,
//                   const double& err, const double& sl_thresh);
// double pdf_dv_log(const double& t, const double& resp, const double& a,
//                   const double& v, const double& w, const double& sv,
//                   const double& err, const double& sl_thresh);
// double pdf_dw_log(const double& t, const double& resp, const double& a,
//                   const double& v, const double& w, const double& sv,
//                   const double& err, const double& sl_thresh);
// double pdf_dsv_log(const double& t, const double& resp, const double& a,
//                    const double& v, const double& w, const double& sv,
//                    const double& err, const double& sl_thresh);



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
                     NumericVector& err,
                     vector<double>& out, const double& rt0);
vector<double> partial_pdf(const ParFunc& parf,
                           const NumericVector& rt,
                           const SEXP& response,
                           const NumericVector& a,
                           const NumericVector& v,
                           const NumericVector& t0,
                           const NumericVector& w,
                           const NumericVector& sv,
                           const NumericVector& sigma,
                           const double& sl_thresh,
                           NumericVector err_tol);
