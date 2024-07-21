// Log-likelihood and gradient function for the Ratcliff DDM

#include "fddm_fit.h"
#include "helper_functions/fitting/declarations.h"



RCPP_MODULE(fddm_fit) {
  Rcpp::class_<fddm_fit>("fddm_fit")
    .constructor<vector<double>, SEXP, vector<MatrixXd>, double, vector<string>, string>("Constructor given response times, responses, model matrices, error tolerance, link functions, and the optimizer")
    .field("rt", &fddm_fit::rt)
    .field("response", &fddm_fit::response)
    .field("err_tol", &fddm_fit::err_tol)
    .field("rt0", &fddm_fit::rt0)
    .field("likelihood", &fddm_fit::likelihood)
    .field("coefficients", &fddm_fit::coefs)
    .field("modmat_v", &fddm_fit::mm_v)
    .field("modmat_a", &fddm_fit::mm_a)
    .field("modmat_t0", &fddm_fit::mm_t0)
    .field("modmat_w", &fddm_fit::mm_w)
    .field("modmat_sv", &fddm_fit::mm_sv)
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
// 5. [optional] list of link functions for the model parameters
// 6. [optional] string name of the optimizer to be used
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
