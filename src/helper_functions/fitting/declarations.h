// Function declarations and constant definitions for fitting the DDM

#ifndef FITTING_DECS_H
#define FITTING_DECS_H

#include <Rcpp.h>
// [[Rcpp::plugins(cpp11)]]
#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]
#include <vector>

#include <cmath>

using std::exp;
using std::log;
using std::sqrt;
using std::max;
using std::isfinite;
using std::isnan;
using std::string;
using std::to_string;
using Rcpp::stop;
using Rcpp::warning;
using std::vector;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::RowVectorXd;
using Eigen::Map;



//--------------------------- Constants --------------------------------------//
static const double ERR_TOL_THRESH = 1e-300; // near minimum value of a double
static const double SV_THRESH = 0; // threshold for using variable drift rate

static const double FOUR_NINTHS = 0.444444444444444444444;
static const double SQRT_2 = 1.41421356237309504880;
static const double SQRT_3 = sqrt(3);
static const double SQRT_5 = sqrt(5);
static const double LOG_2 = 0.693147180559945309417;
static const double LOG_4RT2 = log(4.0 * SQRT_2);
static const double LOG_5_112 = log(5.0/112.0);

static const double PI_CONST = 3.14159265358979323846; // define pi like C++
static const double LOG_PI = log(PI_CONST);
static const double LOG_2PI_2 = 0.5 * log(2 * PI_CONST);
static const double O_PI = 0.318309886183790671538;
static const double PI2 = PI_CONST * PI_CONST;
static const double PI3 = PI_CONST * PI_CONST * PI_CONST;
static const double PI5 = PI_CONST * PI_CONST * PI_CONST * PI_CONST * PI_CONST;
static const double SQRT_2PI = sqrt(2 * PI_CONST);
static const double SQRT_1_2PI = 1 / SQRT_2PI;
static const double SQRT_2_PI = sqrt(2 * O_PI);
static const double SQRT_2_1_PI = SQRT_2 * O_PI;

// static const int INT_MAX = 2147483647; // already in cmath or Rcpp
//----------------------------------------------------------------------------//



//--------------------------- Function Declarations --------------------------//
// Parameter and Input Checks
vector<double> check_rt(const vector<double>& rt, int& Nrt);
vector<double> convert_responses(const SEXP& response, const int& Nrt);
void unpack_and_check_mod_mats(const vector<MatrixXd>& model_matrices,
                               MatrixXd& mm_v, MatrixXd& mm_a, MatrixXd& mm_t0,
                               MatrixXd& mm_w, MatrixXd& mm_sv,
                               VectorXd& v, VectorXd& a, VectorXd& t0,
                               VectorXd& w, VectorXd& sv,
                               vector<int>& form_len, const int& Nrt);
double check_err_tol(const double& err_tol);
double check_switch_thresh(const double& switch_thresh);
bool invalid_parameters(const VectorXd& v, const VectorXd& a,
                        const VectorXd& t0, const VectorXd& w,
                        const VectorXd& sv, const int& Nrt,
                        const vector<int>& form_len);

// PDF (likelihood)
double pdf(const double& t, const double& a, const double& v, const double& w,
           const double& sv, const double& err, const double& switch_thresh);

// Derivatives of PDF
double dv(const double& t, const double& a, const double& v, const double& w,
          const double& sv, const double& err, const double& sl_thresh);
double da(const double& t, const double& a, const double& v, const double& w,
          const double& sv, const double& err, const double& sl_thresh);
double dt(const double& t, const double& a, const double& v, const double& w,
          const double& sv, const double& err, const double& sl_thresh);
double dt0(const double& t, const double& a, const double& v, const double& w,
           const double& sv, const double& err, const double& sl_thresh);
double dw(const double& t, const double& a, const double& v, const double& w,
          const double& sv, const double& err, const double& sl_thresh);
double dsv(const double& t, const double& a, const double& v, const double& w,
           const double& sv, const double& err, const double& sl_thresh);

// Second Order Derivatives of PDF
double dv2(const double& t, const double& a, const double& v, const double& w,
           const double& sv, const double& err, const double& sl_thresh);
double da2(const double& t, const double& a, const double& v, const double& w,
           const double& sv, const double& err, const double& sl_thresh);
double dt02(const double& t, const double& a, const double& v, const double& w,
            const double& sv, const double& err, const double& sl_thresh);
double dw2(const double& t, const double& a, const double& v, const double& w,
           const double& sv, const double& err, const double& sl_thresh);
double dsv2(const double& t, const double& a, const double& v, const double& w,
            const double& sv, const double& err, const double& sl_thresh);

// Infinite Sum Approximations
double small_sum(const double& taa, const double& w, const double& err);
double small_sum_dat(const double& taa, const double& w, const double& err);
double small_sum_dw(const double& taa, const double& w, const int& ks);
double small_sum_dat2(const double& taa, const double& w, const double& err);
double large_sum(const double& taa, const double& w, const int& kl);
double large_sum_dat(const double& taa, const double& w, const int& kl);
double large_sum_dw(const double& taa, const double& w, const int& kl);
double large_sum_dat2(const double& taa, const double& w, const int& kl);

// Required Number of Terms for Infinite Sum Approximations
int kl_pdf(const double& taa, const double& err);
int kl_dat(const double& taa, const double& t, const double& err);
int kl_dw(const double& taa, const double& t, const double& err);
int kl_dat2(const double& taa, const double& err);
int ks_dw(const double& taa, const double& w, const double& err);
//----------------------------------------------------------------------------//

#endif // FITTING_DECS_H
