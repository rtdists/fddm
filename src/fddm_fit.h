// Log-likelihood and gradient function for the Ratcliff DDM

#ifndef FDDM_FIT_H
#define FDDM_FIT_H

#include "helper_functions/fitting/declarations.h"


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
             const double& error_tolerance);

    // methods
    double calc_loglik(const VectorXd& temp_params);
    VectorXd calc_gradient(const VectorXd& temp_params);
    List calc_hessians(const VectorXd& temp_params);
    List calc_vcov();
    VectorXd calc_std_err();
};

#endif // FDDM_FIT_H

// This class requires inputs of:
//
// 1. rt (response times)
// 2. response (associated responses)
// 3. list of model matrices (in order: v, a, t0, w, sv)
// 4. [optional] error_tolerance (for likelihood calculation)
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
