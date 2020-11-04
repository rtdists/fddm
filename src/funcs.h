// Header file to declare constants and functions

#include <Rcpp.h>
// [[Rcpp::plugins(cpp11)]]
#define _USE_MATH_DEFINES
#include <cmath>

using std::vector;
using std::string;
using std::log;
using std::exp;
using std::sqrt;
using std::min;
using std::max;
using std::ceil;
using Rcpp::stop;
using Rcpp::NumericVector;
using Rcpp::LogicalVector;



// Constants

static const double SV_THRESH = 0; // threshold for using variable drift rate
static const double LOG_100 = log(100);
static const double LOG_PI = log(M_PI);
static const double LOG_2PI_2 = 0.5 * log(2 * M_PI);
static const double SQRT_2PI = sqrt(2 * M_PI);
static const char EMPTYCHAR = '\0'; // literally, the empty character
// INT_MAX is in num_funcs.cpp, maximum value of the int data type = 2147483647



// define types for functions in if-else in cpp_dfddm.cpp
typedef int    (*NumFunc)(const double&, const double&, const double&);
typedef double (*SumFunc)(const double&, const double&, const double&,
                          const int&, const double&);
typedef double (*DenFunc)(const double&, const double&, const double&,
                          const double&, const double&, const double&,
                          const int&, const NumFunc&, const SumFunc&);



// Number of Terms

int ks_Gon(const double& t, const double& w, const double& eps);
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
double ff(const double& t, const double& a, const double& v,
          const double& w, const double& sv, const double& eps,
          const int& max_terms_large, const NumFunc& numf, const SumFunc& sumf);
double ff_log(const double& t, const double& a, const double& v,
              const double& w, const double& sv, const double& eps,
              const int& max_terms_large,
              const NumFunc& numf, const SumFunc& sumf);
double fs(const double& t, const double& a, const double& v,
          const double& w, const double& sv, const double& eps,
          const int& max_terms_large, const NumFunc& numf, const SumFunc& sumf);
double fs_log(const double& t, const double& a, const double& v,
              const double& w, const double& sv, const double& eps,
              const int& max_terms_large,
              const NumFunc& numf, const SumFunc& sumf);
double fl(const double& t, const double& a, const double& v,
          const double& w, const double& sv, const double& eps,
          const int& max_terms_large, const NumFunc& numf, const SumFunc& sumf);
double fl_log(const double& t, const double& a, const double& v,
              const double& w, const double& sv, const double& eps,
              const int& max_terms_large,
              const NumFunc& numf, const SumFunc& sumf);
double fc(const double& t, const double& a, const double& v,
          const double& w, const double& sv, const double& eps,
          const int& max_terms_large, const NumFunc& numf, const SumFunc& sumf);
double fc_log(const double& t, const double& a, const double& v,
              const double& w, const double& sv, const double& eps,
              const int& max_terms_large,
              const NumFunc& numf, const SumFunc& sumf);
double fb(const double& t, const double& a, const double& v,
          const double& w, const double& sv, const double& eps,
          const int& max_terms_large, const NumFunc& numf, const SumFunc& sumf);
double fb_log(const double& t, const double& a, const double& v,
              const double& w, const double& sv, const double& eps,
              const int& max_terms_large,
              const NumFunc& numf, const SumFunc& sumf);
