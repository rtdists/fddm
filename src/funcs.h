// Header file to declare constants and functions used in general_density.cpp

#include <Rcpp.h>
// [[Rcpp::plugins(cpp11)]]
#define _USE_MATH_DEFINES
#include <math.h>

// Constants

static const double SV_THRESH = 0.05; // threshold for using variable drift rate
static const double LOG_PI = log(M_PI);
static const double LOG_2PI_2 = 0.5 * log(2 * M_PI);

// Macros (from extraDistr)

#define GETV(x, i) x[i % x.length()] // wrapped indexing of vector

// Number of Terms

int ks_BGK(const double& t, const double& w, const double& eps);
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

Rcpp::NumericVector fs_eps_2017(const Rcpp::NumericVector& rt,
                                Rcpp::LogicalVector response,
                                const double& a, const double& v,
                                const double& t0, const double& w,
                                const double& sv, const double& eps);
Rcpp::NumericVector fs_eps_2014(const Rcpp::NumericVector& rt,
                                Rcpp::LogicalVector response,
                                const double& a, const double& v,
                                const double& t0, const double& w,
                                const double& sv, const double& eps);
Rcpp::NumericVector fs_Nav_2017(const Rcpp::NumericVector& rt,
                                Rcpp::LogicalVector response,
                                const double& a, const double& v,
                                const double& t0, const double& w,
                                const double& sv, const double& eps);
Rcpp::NumericVector fs_Nav_2014(const Rcpp::NumericVector& rt,
                                Rcpp::LogicalVector response,
                                const double& a, const double& v,
                                const double& t0, const double& w,
                                const double& sv, const double& eps);
Rcpp::NumericVector fs_BGK_2017(const Rcpp::NumericVector& rt,
                                Rcpp::LogicalVector response,
                                const double& a, const double& v,
                                const double& t0, const double& w,
                                const double& sv, const double& eps);
Rcpp::NumericVector fs_BGK_2014(const Rcpp::NumericVector& rt,
                                Rcpp::LogicalVector response,
                                const double& a, const double& v,
                                const double& t0, const double& w,
                                const double& sv, const double& eps);
Rcpp::NumericVector fl_Nav_2009(const Rcpp::NumericVector& rt,
                                Rcpp::LogicalVector response,
                                const double& a, const double& v,
                                const double& t0, const double& w,
                                const double& eps);
Rcpp::NumericVector fb_Nav_Nav_2017(const Rcpp::NumericVector& rt,
                                    Rcpp::LogicalVector response,
                                    const double& a, const double& v,
                                    const double& t0, const double& w,
                                    const double& sv, const double& eps);
Rcpp::NumericVector fb_Nav_Nav_2014(const Rcpp::NumericVector& rt,
                                    Rcpp::LogicalVector response,
                                    const double& a, const double& v,
                                    const double& t0, const double& w,
                                    const double& sv, const double& eps);
Rcpp::NumericVector fb_BGK_Nav_2017(const Rcpp::NumericVector& rt,
                                    Rcpp::LogicalVector response,
                                    const double& a, const double& v,
                                    const double& t0, const double& w,
                                    const double& sv, const double& eps);
Rcpp::NumericVector fb_BGK_Nav_2014(const Rcpp::NumericVector& rt,
                                    Rcpp::LogicalVector response,
                                    const double& a, const double& v,
                                    const double& t0, const double& w,
                                    const double& sv, const double& eps);

// Density Functions (log)

Rcpp::NumericVector fs_eps_2017_log(const Rcpp::NumericVector& rt,
                                    Rcpp::LogicalVector response,
                                    const double& a, const double& v,
                                    const double& t0, const double& w,
                                    const double& sv, const double& eps);
Rcpp::NumericVector fs_eps_2014_log(const Rcpp::NumericVector& rt,
                                    Rcpp::LogicalVector response,
                                    const double& a, const double& v,
                                    const double& t0, const double& w,
                                    const double& sv, const double& eps);
Rcpp::NumericVector fs_Nav_2017_log(const Rcpp::NumericVector& rt,
                                    Rcpp::LogicalVector response,
                                    const double& a, const double& v,
                                    const double& t0, const double& w,
                                    const double& sv, const double& eps);
Rcpp::NumericVector fs_Nav_2014_log(const Rcpp::NumericVector& rt,
                                    Rcpp::LogicalVector response,
                                    const double& a, const double& v,
                                    const double& t0, const double& w,
                                    const double& sv, const double& eps);
Rcpp::NumericVector fs_BGK_2017_log(const Rcpp::NumericVector& rt,
                                    Rcpp::LogicalVector response,
                                    const double& a, const double& v,
                                    const double& t0, const double& w,
                                    const double& sv, const double& eps);
Rcpp::NumericVector fs_BGK_2014_log(const Rcpp::NumericVector& rt,
                                    Rcpp::LogicalVector response,
                                    const double& a, const double& v,
                                    const double& t0, const double& w,
                                    const double& sv, const double& eps);
Rcpp::NumericVector fl_Nav_2009_log(const Rcpp::NumericVector& rt,
                                    Rcpp::LogicalVector response,
                                    const double& a, const double& v,
                                    const double& t0, const double& w,
                                    const double& eps);
Rcpp::NumericVector fb_Nav_Nav_2017_log(const Rcpp::NumericVector& rt,
                                        Rcpp::LogicalVector response,
                                        const double& a, const double& v,
                                        const double& t0, const double& w,
                                        const double& sv, const double& eps);
Rcpp::NumericVector fb_Nav_Nav_2014_log(const Rcpp::NumericVector& rt,
                                        Rcpp::LogicalVector response,
                                        const double& a, const double& v,
                                        const double& t0, const double& w,
                                        const double& sv, const double& eps);
Rcpp::NumericVector fb_BGK_Nav_2017_log(const Rcpp::NumericVector& rt,
                                        Rcpp::LogicalVector response,
                                        const double& a, const double& v,
                                        const double& t0, const double& w,
                                        const double& sv, const double& eps);
Rcpp::NumericVector fb_BGK_Nav_2014_log(const Rcpp::NumericVector& rt,
                                        Rcpp::LogicalVector response,
                                        const double& a, const double& v,
                                        const double& t0, const double& w,
                                        const double& sv, const double& eps);
