// Header file to declare constants and functions used in general_density.cpp

#include <Rcpp.h>
// [[Rcpp::plugins(cpp11)]]
#define _USE_MATH_DEFINES
#include <cmath>

using std::log;
using std::exp;
using std::sqrt;
using std::min;
using std::max;
using std::ceil;
using Rcpp::NumericVector;
using Rcpp::LogicalVector;



// Constants

static const double SV_THRESH = 0.05; // threshold for using variable drift rate
static const double LOG_PI = log(M_PI);
static const double LOG_2PI_2 = 0.5 * log(2 * M_PI);



// define type for density function in if-else in cpp_dfddm.cpp
typedef double (*DensFunc)(const double&, const double&, const double&,
                           const double&, const double&, const bool&,
                           const double&);



// Number of Terms

int ks_Kes(const double& t, const double& w, const double& eps);
int ks_Nav(const double& t, const double& eps);
int kl_Nav(const double& t, const double& eps);



// Infinite Summation Approximations

double small_sum_eps_17(const double& t, const double& a, const double& w,
                        const double& eps);
double small_sum_eps_14(const double& t, const double& a, const double& w,
                        const double&  eps);
double small_sum_2017(const double& t, const double& a, const double& w,
                      const int& ks);
double small_sum_2014(const double& t, const double& a, const double& w,
                      const int& ks);
double large_sum_Nav(const double& t, const double& a, const double& w,
                     const int& kl);



// Density Functions

double fs_Fos_17(const double& t, const double& a,  const double& v,
                 const double& w, const double& sv, const bool& log_prob,
                 const double& eps);
double fs_Fos_14(const double& t, const double& a,  const double& v,
                 const double& w, const double& sv, const bool& log_prob,
                 const double& eps);
double fs_Kes_17(const double& t, const double& a,  const double& v,
                 const double& w, const double& sv, const bool& log_prob,
                 const double& eps);
double fs_Kes_14(const double& t, const double& a,  const double& v,
                 const double& w, const double& sv, const bool& log_prob,
                 const double& eps);
double fs_Nav_17(const double& t, const double& a,  const double& v,
                 const double& w, const double& sv, const bool& log_prob,
                 const double& eps);
double fs_Nav_14(const double& t, const double& a,  const double& v,
                 const double& w, const double& sv, const bool& log_prob,
                 const double& eps);
double fl_Nav_09(const double& t, const double& a,  const double& v,
                 const double& w, const double& sv, const bool& log_prob,
                 const double& eps);
double fb_Kes_17(const double& t, const double& a,  const double& v,
                 const double& w, const double& sv, const bool& log_prob,
                 const double& eps);
double fb_Kes_14(const double& t, const double& a,  const double& v,
                 const double& w, const double& sv, const bool& log_prob,
                 const double& eps);
double fb_Nav_17(const double& t, const double& a,  const double& v,
                 const double& w, const double& sv, const bool& log_prob,
                 const double& eps);
double fb_Nav_14(const double& t, const double& a,  const double& v,
                 const double& w, const double& sv, const bool& log_prob,
                 const double& eps);
