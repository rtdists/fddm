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

static const double SV_THRESH = 1e-6; // threshold for using variable drift rate
static const double LOG_PI = log(M_PI);
static const double LOG_2PI_2 = 0.5 * log(2 * M_PI);
static const double SQRT_2PI = sqrt(2 * M_PI);



// define types for functions in if-else in cpp_dfddm.cpp
typedef int    (*NummFunc)(const double&, const double&, const double&);
typedef double (*SummFunc)(const double&, const double&, const double&,
                           const int&, const double&);
typedef double (*DensFunc)(const double&, const double&, const double&,
                           const double&, const double&, const bool&,
                           const double&, NummFunc, SummFunc, int);



// Number of Terms

int ks_Kes(const double& t, const double& w, const double& eps);
int ks_Nav(const double& t, const double& w, const double& eps);
int kl_Nav(const double& t, const double& w, const double& eps);



// Infinite Summation Approximations

double small_sum_eps_17(const double& t, const double& a, const double& w,
                        const int& ks, const double& eps);
double small_sum_eps_14(const double& t, const double& a, const double& w,
                        const int& ks, const double&  eps);
double small_sum_2017(const double& t, const double& a, const double& w,
                      const int& ks, const double& eps);
double small_sum_2014(const double& t, const double& a, const double& w,
                      const int& ks, const double& eps);
double large_sum_Nav(const double& t, const double& a, const double& w,
                     const int& kl, const double& eps);



// Density Functions
double fs(const double& t, const double& a, const double& v, const double& w,
          const double& sv, const bool& log_prob, const double& eps,
          NummFunc numm, SummFunc summ, int ks);
double fl(const double& t, const double& a, const double& v, const double& w,
          const double& sv, const bool& log_prob, const double& eps,
          NummFunc numm, SummFunc summ, int kl);
double fb(const double& t, const double& a, const double& v, const double& w,
          const double& sv, const bool& log_prob, const double& eps,
          NummFunc numm, SummFunc summ, int k);
