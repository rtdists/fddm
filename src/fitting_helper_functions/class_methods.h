// included from src/fitting_class.cpp

#include "declarations.h" // all function declarations and constant definitions

#include "checks.h" // src/fitting_helper_functions/par_checks.h
#include "likelihood_funcs.h" // src/fitting_helper_functions/likelihood_funcs.h
#include "gradient_funcs.h" // src/fitting_helper_functions/gradient_funcs.h
#include "sum_funcs.h" // src/fitting_helper_functions/sum_funcs.h
#include "num_funcs.h" // src/fitting_helper_functions/num_funcs.h

using std::log;
using std::isfinite;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::Map;



// Constructors
fddm_fit::fddm_fit(const vector<double>& rt_vector,
                   const SEXP& response_vector,
                   const vector<MatrixXd>& model_matrices,
                   const double& error_tolerance,
                   const double& switching_threshold)
  {
  rt = check_rt(rt_vector, Nrt); // also sets Nrt
  response = convert_responses(response_vector, Nrt);
  likelihood.resize(Nrt);
  v.resize(Nrt);
  a.resize(Nrt);
  t0.resize(Nrt);
  w.resize(Nrt);
  sv.resize(Nrt);
  unpack_and_check_mod_mats(model_matrices, mm_v, mm_a, mm_t0, mm_w, mm_sv,
                            v, a, t0, w, sv, form_len, Nrt);
  err_tol = check_err_tol(error_tolerance);
  switch_thresh = check_switch_thresh(switching_threshold);
  // Resize coefs
  for (int i = 0; i < 5; i++) { Ncoefs += form_len[i]; }
  coefs.resize(Ncoefs);
}


// Log-likelihood Function
double fddm_fit::calc_loglik(const VectorXd& temp_coefs)
{
  // Store the input parameters for later checking (in gradient)
  coefs = temp_coefs;

  // Split input parameter vector and multiply with model matrix
  int coef_idx = 0;
  if (form_len[0] > 0) { // check if v has a model matrix (i.e., is a formula)
    v = mm_v * coefs.segment(coef_idx, form_len[0]);
    coef_idx += form_len[0];
  } // else it's constant and was handled in the constructor

  if (form_len[1] > 0) { // check if a has a model matrix (i.e., is a formula)
    a = mm_a * coefs.segment(coef_idx, form_len[1]);
    coef_idx += form_len[1];
  } // else it's constant and was handled in the constructor

  if (form_len[2] > 0) { // check if t0 has a model matrix (i.e., is a formula)
    t0 = mm_t0 * coefs.segment(coef_idx, form_len[2]);
    coef_idx += form_len[2];
  } // else it's constant and was handled in the constructor

  if (form_len[3] > 0) { // check if w has a model matrix (i.e., is a formula)
    w = mm_w * coefs.segment(coef_idx, form_len[3]);
    coef_idx += form_len[3];
  } // else it's constant and was handled in the constructor

  if (form_len[4] > 0) { // check if sv has a model matrix (i.e., is a formula)
    sv = mm_sv * coefs.segment(coef_idx, form_len[4]);
  } // else it's constant and was handled in the constructor

  // Check parameters
  if (invalid_parameters(v, a, t0, w, sv, Nrt, form_len)) {
    for (int i = 0; i < Nrt; i++) {
      likelihood[i] = rt0;
    }
    return rt0;
  }

  // Calculate loglikelihood and -sum over all data
  double ll = 0.0;
  double t;
  for (int i = 0; i < Nrt; i++) {
    t = rt[i] - t0[i];
    if (t > 0 && isfinite(t)) {
      if (response[i] == 1) { // response is "lower"
        likelihood[i] = ft(t, a[i], v[i], w[i], sv[i], err_tol, switch_thresh);
      } else { // response is "upper" so use alternate parameters
        likelihood[i] = ft(t, a[i], -v[i], 1-w[i], sv[i], err_tol,
                           switch_thresh);
      }
      ll -= log(likelihood[i]); // faster than doing exp(loglik)
    } else {
      likelihood[i] = rt0;
      return rt0;
    }
  }
  return ll;
}

// Log-gradient Function
VectorXd fddm_fit::calc_gradient(const VectorXd& temp_coefs)
{
  // check if new coefficients match the stored ones
  if (temp_coefs != coefs) { // if not, must recalculate likelihood
    double not_used = calc_loglik(temp_coefs); // stored in `likelihood`
  } // if they match, can reuse the likelihood
  
  // calculate gradient and -sum over all data
  // note that d/dx(log(f(x))) = d/dx(f(x)) / f(x) (why dividing by likelihood)
  VectorXd gradient = VectorXd::Zero(Ncoefs);
  int grad_idx;
  double t, temp_grad;
  for (int i = 0; i < Nrt; i++) {
    grad_idx = 0;
    t = rt[i] - t0[i];
    if (t > 0 && isfinite(t)) { // maybe don't have to check b/c same in lik calc?
      if (response[i] == 1) { // response is "lower"
        if (form_len[0] > 0) {
          temp_grad = dv(t, a[i], v[i], w[i], sv[i], err_tol, switch_thresh)
                      / likelihood[i];
          gradient.segment(grad_idx, form_len[0]) -= temp_grad * mm_v.row(i);
          grad_idx += form_len[0];
        }
        if (form_len[1] > 0) {
          temp_grad = da(t, a[i], v[i], w[i], sv[i], err_tol, switch_thresh)
                      / likelihood[i];
          gradient.segment(grad_idx, form_len[1]) -= temp_grad * mm_a.row(i);
          grad_idx += form_len[1];
        }
        if (form_len[2] > 0) {
          temp_grad = dt0(t, a[i], v[i], w[i], sv[i], err_tol, switch_thresh)
                      / likelihood[i];
          gradient.segment(grad_idx, form_len[2]) -= temp_grad * mm_t0.row(i);
          grad_idx += form_len[2];
        }
        if (form_len[3] > 0) {
          temp_grad = dw(t, a[i], v[i], w[i], sv[i], err_tol, switch_thresh)
                      / likelihood[i];
          gradient.segment(grad_idx, form_len[3]) -= temp_grad * mm_w.row(i);
          grad_idx += form_len[3];
        }
        if (form_len[4] > 0) {
          temp_grad = dsv(t, a[i], v[i], w[i], sv[i], err_tol, switch_thresh)
                      / likelihood[i];
          gradient.segment(grad_idx, form_len[4]) -= temp_grad * mm_sv.row(i);
          grad_idx += form_len[4];
        }
      } else { // response is "upper" so use alternate parameters
        if (form_len[0] > 0) { // chain rule negates derivative
          temp_grad = dv(t, a[i], -v[i], 1-w[i], sv[i], err_tol, switch_thresh)
                      / likelihood[i];
          gradient.segment(grad_idx, form_len[0]) += temp_grad * mm_v.row(i);
          grad_idx += form_len[0];
        }
        if (form_len[1] > 0) {
          temp_grad = da(t, a[i], -v[i], 1-w[i], sv[i], err_tol, switch_thresh)
                      / likelihood[i];
          gradient.segment(grad_idx, form_len[1]) -= temp_grad * mm_a.row(i);
          grad_idx += form_len[1];
        }
        if (form_len[2] > 0) {
          temp_grad = dt0(t, a[i], -v[i], 1-w[i], sv[i], err_tol, switch_thresh)
                      / likelihood[i];
          gradient.segment(grad_idx, form_len[2]) -= temp_grad * mm_t0.row(i);
          grad_idx += form_len[2];
        }
        if (form_len[3] > 0) { // chain rule negates derivative
          temp_grad = dw(t, a[i], -v[i], 1-w[i], sv[i], err_tol, switch_thresh)
                      / likelihood[i];
          gradient.segment(grad_idx, form_len[3]) += temp_grad * mm_w.row(i);
          grad_idx += form_len[3];
        }
        if (form_len[4] > 0) {
          temp_grad = dsv(t, a[i], -v[i], 1-w[i], sv[i], err_tol, switch_thresh)
                      / likelihood[i];
          gradient.segment(grad_idx, form_len[4]) -= temp_grad * mm_sv.row(i);
          grad_idx += form_len[4];
        }
      }
    } else {
      likelihood[i] = rt0;
      for (int j = 0; j < Ncoefs; j++) {
        gradient[j] = std::numeric_limits<double>::quiet_NaN();
      }
      break;
    }
  }
  return gradient;
}
