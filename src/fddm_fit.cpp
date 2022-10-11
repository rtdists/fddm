// Log-likelihood and gradient function for the Ratcliff DDM
#include <Rcpp.h>
// [[Rcpp::plugins(cpp11)]]
#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]
#include <vector>

using std::vector;
using Eigen::MatrixXd;
using Eigen::VectorXd;

//' @export fddm_fit
class fddm_fit {
  public:
    // variables (accessible from R)
    vector<double> rt {};
    vector<double> response {};
    MatrixXd mm_v {};
    MatrixXd mm_a {};
    MatrixXd mm_t0 {};
    MatrixXd mm_w {};
    MatrixXd mm_sv {};
    double err_tol {0.000001};
    double switch_thresh {0.8};
    vector<double> likelihood {};
    VectorXd coefs {};
    MatrixXd hess_v {};
    MatrixXd hess_a {};
    MatrixXd hess_t0 {};
    MatrixXd hess_w {};
    MatrixXd hess_sv {};
    MatrixXd vcov_v {};
    MatrixXd vcov_a {};
    MatrixXd vcov_t0 {};
    MatrixXd vcov_w {};
    MatrixXd vcov_sv {};

    // variables (internal use only)
    int Nrt {};
    double rt0 {1e6};
    vector<int> form_len {0, 0, 0, 0, 0}; // length of each parameter's formula
    int Ncoefs {0};
    VectorXd v {};
    VectorXd a {};
    VectorXd t0 {};
    VectorXd w {};
    VectorXd sv {};

    // constructors
    fddm_fit(const vector<double>& rt_vector,
             const SEXP& response_vector,
             const vector<MatrixXd>& model_matrices,
             const double& error_tolerance,
             const double& switching_threshold);

    // methods
    double calc_loglik(const VectorXd& temp_params);
    VectorXd calc_gradient(const VectorXd& temp_params);
    void calc_hessians(const VectorXd& temp_params);
    void calc_vcov();
    VectorXd calc_std_err();
};

#include "fitting_helper_functions/class_methods.h"


RCPP_MODULE(fddm_fit) {
  using namespace Rcpp;

  class_<fddm_fit>( "fddm_fit")
    .constructor<vector<double>, SEXP, vector<MatrixXd>, double, double>("Constructor given response times, responses, model matrices, error tolerance, and switching threshold")
    .field("rt", &fddm_fit::rt)
    .field("response", &fddm_fit::response)
    .field("modmat_v", &fddm_fit::mm_v)
    .field("modmat_a", &fddm_fit::mm_a)
    .field("modmat_t0", &fddm_fit::mm_t0)
    .field("modmat_w", &fddm_fit::mm_w)
    .field("modmat_sv", &fddm_fit::mm_sv)
    .field("err_tol", &fddm_fit::err_tol)
    .field("switch_thresh", &fddm_fit::switch_thresh)
    .field("likelihood", &fddm_fit::likelihood)
    .field("coefficients", &fddm_fit::coefs)
    .field("hess_v", &fddm_fit::hess_v)
    .field("hess_a", &fddm_fit::hess_a)
    .field("hess_t0", &fddm_fit::hess_t0)
    .field("hess_w", &fddm_fit::hess_w)
    .field("hess_sv", &fddm_fit::hess_sv)
    .field("vcov_v", &fddm_fit::vcov_v)
    .field("vcov_a", &fddm_fit::vcov_a)
    .field("vcov_t0", &fddm_fit::vcov_t0)
    .field("vcov_w", &fddm_fit::vcov_w)
    .field("vcov_sv", &fddm_fit::vcov_sv)
    .method("calculate_loglik", &fddm_fit::calc_loglik)
    .method("calculate_gradient", &fddm_fit::calc_gradient)
    .method("calculate_hessians", &fddm_fit::calc_hessians)
    .method("calculate_vcov", &fddm_fit::calc_vcov)
    .method("calculate_standard_error", &fddm_fit::calc_std_err)
  ;
}



// This class requires inputs of:
//
// 1. rt (response times)
// 2. response (associated responses)
// 3. list of model matrices (in order: v, a, t0, w, sv)
// 4. [optional] error_tolerance (for likelihood calculation)
// 5. [optional] switching threshold (for likelihood calculation)
//
//
// Note on model matrices:
// The model matrix can be a 1x1 matrix with a single value. In this case, the
// respective DDM parameter is treated as a constant, and the single value in
// the 1x1 model matrix will be used as a constant value during the fitting.
// The model matrices can also be an nx1 column vector, indicating a formula
// with a single feature (i.e., ~ 1 or ~ 0 + feat). Note that if the data has
// only one observation, this will result in the model matrix being read as a
// constant.
