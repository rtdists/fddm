// Header file to declare constants and functions used in general_density.cpp

#include <Rcpp.h>
// [[Rcpp::plugins(cpp11)]]
#define _USE_MATH_DEFINES
#include <math.h>

// Constants

static const double SV_THRESH = 0.05; // threshold for using variable drift rate

// Number of Terms

int ks_BGK(double t, double a, double w, double eps);
int ks_Nav(double t, double a, double eps);
int kl_Nav(double t, double a, double eps);

// Infinite Summation Approximations

double small_sum_2014(double t, double a, double w, int ks);
double small_sum_2017(double t, double a, double w, int ks);
double small_sum_eps_14(double t, double a, double w, double mult, double eps);
double small_sum_eps_17(double t, double a, double w, double mult, double eps);
double large_sum_Nav(double t, double a, double w, int kl);

// Density Functions

Rcpp::NumericVector fs_eps_2017(Rcpp::NumericVector rt,
                                Rcpp::LogicalVector response,
                                double a, double v, double t0, double w,
                                double sv, bool log, double eps);
Rcpp::NumericVector fs_eps_2014(Rcpp::NumericVector rt,
                                Rcpp::LogicalVector response,
                                double a, double v, double t0, double w,
                                double sv, bool log, double eps);
Rcpp::NumericVector fs_Nav_2017(Rcpp::NumericVector rt,
                                Rcpp::LogicalVector response,
                                double a, double v, double t0, double w,
                                double sv, bool log, double eps);
Rcpp::NumericVector fs_Nav_2014(Rcpp::NumericVector rt,
                                Rcpp::LogicalVector response,
                                double a, double v, double t0, double w,
                                double sv, bool log, double eps);
Rcpp::NumericVector fs_BGK_2017(Rcpp::NumericVector rt,
                                Rcpp::LogicalVector response,
                                double a, double v, double t0, double w,
                                double sv, bool log, double eps);
Rcpp::NumericVector fs_BGK_2014(Rcpp::NumericVector rt,
                                Rcpp::LogicalVector response,
                                double a, double v, double t0, double w,
                                double sv, bool log, double eps);
Rcpp::NumericVector fl_Nav_2009(Rcpp::NumericVector rt,
                                Rcpp::LogicalVector response,
                                double a, double v, double t0, double w,
                                bool log, double eps);
Rcpp::NumericVector fb_Nav_Nav_2017(Rcpp::NumericVector rt,
                                    Rcpp::LogicalVector response,
                                    double a, double v, double t0, double w,
                                    double sv, bool log_prob, double eps);
Rcpp::NumericVector fb_Nav_Nav_2014(Rcpp::NumericVector rt,
                                    Rcpp::LogicalVector response,
                                    double a, double v, double t0, double w,
                                    double sv, bool log_prob, double eps);
Rcpp::NumericVector fb_BGK_Nav_2017(Rcpp::NumericVector rt,
                                    Rcpp::LogicalVector response,
                                    double a, double v, double t0, double w,
                                    double sv, bool log_prob, double eps);
Rcpp::NumericVector fb_BGK_Nav_2014(Rcpp::NumericVector rt,
                                    Rcpp::LogicalVector response,
                                    double a, double v, double t0, double w,
                                    double sv, bool log_prob, double eps);
