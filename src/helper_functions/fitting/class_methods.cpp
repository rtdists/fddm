// Methods for the fddm_fit class

#include "../../fddm_fit.h"
#include "declarations.h" // namespaces, function declarations, constant definitions



// Constructors
fddm_fit::fddm_fit(const vector<double>& rt_vector,
                   const SEXP& response_vector,
                   const vector<MatrixXd>& model_matrices,
                   const double& error_tolerance)
  {
  rt = check_rt(rt_vector, Nrt, min_rt); // also sets Nrt and min_rt
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
  // Resize coefs
  for (int i = 0; i < 5; i++) { Ncoefs += form_len[i]; }
  coefs.resize(Ncoefs);
}


// Log-likelihood Function (negated)
double fddm_fit::calc_loglik(const VectorXd& temp_coefs)
{
  // Reset likelihood flag (for invalid parameters)
  lik_flag = 0;
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
  if (invalid_parameters(v, a, t0, w, sv, Nrt, min_rt, form_len)) {
    for (int i = 0; i < Nrt; i++) {
      likelihood[i] = rt0;
    }
    lik_flag = 1;
    return rt0;
  }

  // Calculate (negative) loglikelihood and (-)sum over all data
  double ll = 0.0;
  double t;
  for (int i = 0; i < Nrt; i++) {
    t = rt[i] - t0[i];
    if (t > 0 && isfinite(t)) {
      if (response[i] == 1) { // response is "lower"
        likelihood[i] = pdf(t, v[i], a[i], w[i], sv[i], err_tol);
      } else { // response is "upper" so use alternate parameters
        likelihood[i] = pdf(t, -v[i], a[i], 1-w[i], sv[i], err_tol);
      }
      ll -= log(likelihood[i]); // faster than doing exp(loglik)
    } else {
      likelihood[i] = rt0;
      lik_flag = 1;
      return rt0;
    }
  }
  return ll;
}


// Log-gradient Function (negated)
VectorXd fddm_fit::calc_gradient(const VectorXd& temp_coefs)
{
  // check if new coefficients match the stored ones
  if (temp_coefs != coefs) { // if not, must recalculate likelihood
    double not_used = calc_loglik(temp_coefs); // stored in `likelihood`
  } // if they match, can reuse the likelihood

  // calculate (negative) gradient and (-)sum over all data
  // note that d/dx(log(f(x))) = d/dx(f(x)) / f(x) (why dividing by likelihood)
  VectorXd gradient = VectorXd::Zero(Ncoefs);
  if (lik_flag) { // check if likelihood is rt0 (i.e., a model parameter is invalid)
    return gradient; // or rt0 something like that
  }
  int grad_idx;
  double t, temp_grad;
  for (int i = 0; i < Nrt; i++) {
    grad_idx = 0;
    t = rt[i] - t0[i];
    if (t > 0 && isfinite(t)) { // maybe don't have to check b/c same in lik calc?
      if (response[i] == 1) { // response is "lower"
        if (form_len[0] > 0) {
          temp_grad = dv(t, v[i], a[i], w[i], sv[i], err_tol)
                      / likelihood[i];
          gradient.segment(grad_idx, form_len[0]) -= temp_grad * mm_v.row(i);
          grad_idx += form_len[0];
        }
        if (form_len[1] > 0) {
          temp_grad = da(t, v[i], a[i], w[i], sv[i], err_tol)
                      / likelihood[i];
          gradient.segment(grad_idx, form_len[1]) -= temp_grad * mm_a.row(i);
          grad_idx += form_len[1];
        }
        if (form_len[2] > 0) {
          temp_grad = dt0(t, v[i], a[i], w[i], sv[i], err_tol)
                      / likelihood[i];
          gradient.segment(grad_idx, form_len[2]) -= temp_grad * mm_t0.row(i);
          grad_idx += form_len[2];
        }
        if (form_len[3] > 0) {
          temp_grad = dw(t, v[i], a[i], w[i], sv[i], err_tol)
                      / likelihood[i];
          gradient.segment(grad_idx, form_len[3]) -= temp_grad * mm_w.row(i);
          grad_idx += form_len[3];
        }
        if (form_len[4] > 0) {
          temp_grad = dsv(t, v[i], a[i], w[i], sv[i], err_tol)
                      / likelihood[i];
          gradient.segment(grad_idx, form_len[4]) -= temp_grad * mm_sv.row(i);
          grad_idx += form_len[4];
        }
      } else { // response is "upper" so use alternate parameters
        if (form_len[0] > 0) { // chain rule negates derivative
          temp_grad = dv(t, -v[i], a[i], 1-w[i], sv[i], err_tol)
                      / likelihood[i];
          gradient.segment(grad_idx, form_len[0]) += temp_grad * mm_v.row(i);
          grad_idx += form_len[0];
        }
        if (form_len[1] > 0) {
          temp_grad = da(t, -v[i], a[i], 1-w[i], sv[i], err_tol)
                      / likelihood[i];
          gradient.segment(grad_idx, form_len[1]) -= temp_grad * mm_a.row(i);
          grad_idx += form_len[1];
        }
        if (form_len[2] > 0) {
          temp_grad = dt0(t, -v[i], a[i], 1-w[i], sv[i], err_tol)
                      / likelihood[i];
          gradient.segment(grad_idx, form_len[2]) -= temp_grad * mm_t0.row(i);
          grad_idx += form_len[2];
        }
        if (form_len[3] > 0) { // chain rule negates derivative
          temp_grad = dw(t, -v[i], a[i], 1-w[i], sv[i], err_tol)
                      / likelihood[i];
          gradient.segment(grad_idx, form_len[3]) += temp_grad * mm_w.row(i);
          grad_idx += form_len[3];
        }
        if (form_len[4] > 0) {
          temp_grad = dsv(t, -v[i], a[i], 1-w[i], sv[i], err_tol)
                      / likelihood[i];
          gradient.segment(grad_idx, form_len[4]) -= temp_grad * mm_sv.row(i);
          grad_idx += form_len[4];
        }
      }
    } else {
      likelihood[i] = rt0;
      for (int j = 0; j < Ncoefs; j++) {
        // gradient[j] = std::numeric_limits<double>::quiet_NaN();
        gradient[j] = 0.0; // 0? rt0? I don't know what to put here
      }
      break;
    }
  }
  return gradient;
}


// Log-Hessian Function (negated)
List fddm_fit::calc_hessians(const VectorXd& temp_coefs)
{
  // if MLE coefficients do not match the stored ones, must recalculate
  // parameter vectors, likelihood, gradient
  if (temp_coefs != coefs) {
    double not_used = calc_loglik(temp_coefs); // stored in `likelihood`
  }
  
  // calculate (negative) hessians and (-)sum over all data
  // the second derivative of log(f(x)) is a mess, see paper for details
  hess_v = MatrixXd::Zero(form_len[0], form_len[0]);
  hess_a = MatrixXd::Zero(form_len[1], form_len[1]);
  hess_t0 = MatrixXd::Zero(form_len[2], form_len[2]);
  hess_w = MatrixXd::Zero(form_len[3], form_len[3]);
  hess_sv = MatrixXd::Zero(form_len[4], form_len[4]);
  RowVectorXd mm_v_row = RowVectorXd::Zero(form_len[0]);
  VectorXd mm_v_col = VectorXd::Zero(form_len[0]);
  RowVectorXd mm_a_row = RowVectorXd::Zero(form_len[1]);
  VectorXd mm_a_col = VectorXd::Zero(form_len[1]);
  RowVectorXd mm_t0_row = RowVectorXd::Zero(form_len[2]);
  VectorXd mm_t0_col = VectorXd::Zero(form_len[2]);
  RowVectorXd mm_w_row = RowVectorXd::Zero(form_len[3]);
  VectorXd mm_w_col = VectorXd::Zero(form_len[3]);
  RowVectorXd mm_sv_row = RowVectorXd::Zero(form_len[4]);
  VectorXd mm_sv_col = VectorXd::Zero(form_len[4]);
  double t, d1, d2, dd;
  for (int i = 0; i < Nrt; i++) {
    t = rt[i] - t0[i];
    if (t > 0 && isfinite(t)) { // maybe don't have to check b/c same in lik calc?
      if (response[i] == 1) { // response is "lower"
        if (form_len[0] > 0) {
          d1 = dv(t, v[i], a[i], w[i], sv[i], err_tol);
          d2 = dv2(t, v[i], a[i], w[i], sv[i], err_tol);
          dd = (likelihood[i] * d2 - d1*d1) / (likelihood[i]*likelihood[i]);
          mm_v_row = mm_v.row(i);
          mm_v_col = mm_v.row(i);
          hess_v -= mm_v_col * dd * mm_v_row;
        }
        if (form_len[1] > 0) {
          d1 = da(t, v[i], a[i], w[i], sv[i], err_tol);
          d2 = da2(t, v[i], a[i], w[i], sv[i], err_tol);
          dd = (likelihood[i] * d2 - d1*d1) / (likelihood[i]*likelihood[i]);
          mm_a_row = mm_a.row(i);
          mm_a_col = mm_a.row(i);
          hess_a -= mm_a_col * dd * mm_a_row;
        }
        if (form_len[2] > 0) {
          d1 = dt0(t, v[i], a[i], w[i], sv[i], err_tol);
          d2 = dt02(t, v[i], a[i], w[i], sv[i], err_tol);
          dd = (likelihood[i] * d2 - d1*d1) / (likelihood[i]*likelihood[i]);
          mm_t0_row = mm_t0.row(i);
          mm_t0_col = mm_t0.row(i);
          hess_t0 -= mm_t0_col * dd * mm_t0_row;
        }
        if (form_len[3] > 0) {
          d1 = dw(t, v[i], a[i], w[i], sv[i], err_tol);
          d2 = dw2(t, v[i], a[i], w[i], sv[i], err_tol);
          dd = (likelihood[i] * d2 - d1*d1) / (likelihood[i]*likelihood[i]);
          mm_w_row = mm_w.row(i);
          mm_w_col = mm_w.row(i);
          hess_w -= mm_w_col * dd * mm_w_row;
        }
        if (form_len[4] > 0) {
          d1 = dsv(t, v[i], a[i], w[i], sv[i], err_tol);
          d2 = dsv2(t, v[i], a[i], w[i], sv[i], err_tol);
          dd = (likelihood[i] * d2 - d1*d1) / (likelihood[i]*likelihood[i]);
          mm_sv_row = mm_sv.row(i);
          mm_sv_col = mm_sv.row(i);
          hess_sv -= mm_sv_col * dd * mm_sv_row;
        }
      } else { // response is "upper" so use alternate parameters
        // note: chain rule negation not needed because squared or negated twice
        if (form_len[0] > 0) {
          d1 = dv(t, -v[i], a[i], 1-w[i], sv[i], err_tol);
          d2 = dv2(t, -v[i], a[i], 1-w[i], sv[i], err_tol);
          dd = (likelihood[i] * d2 - d1*d1) / (likelihood[i]*likelihood[i]);
          mm_v_row = mm_v.row(i);
          mm_v_col = mm_v.row(i);
          hess_v -= mm_v_col * dd * mm_v_row;
        }
        if (form_len[1] > 0) {
          d1 = da(t, -v[i], a[i], 1-w[i], sv[i], err_tol);
          d2 = da2(t, -v[i], a[i], 1-w[i], sv[i], err_tol);
          dd = (likelihood[i] * d2 - d1*d1) / (likelihood[i]*likelihood[i]);
          mm_a_row = mm_a.row(i);
          mm_a_col = mm_a.row(i);
          hess_a -= mm_a_col * dd * mm_a_row;
        }
        if (form_len[2] > 0) {
          d1 = dt0(t, -v[i], a[i], 1-w[i], sv[i], err_tol);
          d2 = dt02(t, -v[i], a[i], 1-w[i], sv[i], err_tol);
          dd = (likelihood[i] * d2 - d1*d1) / (likelihood[i]*likelihood[i]);
          mm_t0_row = mm_t0.row(i);
          mm_t0_col = mm_t0.row(i);
          hess_t0 -= mm_t0_col * dd * mm_t0_row;
        }
        if (form_len[3] > 0) {
          d1 = dw(t, -v[i], a[i], 1-w[i], sv[i], err_tol);
          d2 = dw2(t, -v[i], a[i], 1-w[i], sv[i], err_tol);
          dd = (likelihood[i] * d2 - d1*d1) / (likelihood[i]*likelihood[i]);
          mm_w_row = mm_w.row(i);
          mm_w_col = mm_w.row(i);
          hess_w -= mm_w_col * dd * mm_w_row;
        }
        if (form_len[4] > 0) {
          d1 = dsv(t, -v[i], a[i], 1-w[i], sv[i], err_tol);
          d2 = dsv2(t, -v[i], a[i], 1-w[i], sv[i], err_tol);
          dd = (likelihood[i] * d2 - d1*d1) / (likelihood[i]*likelihood[i]);
          mm_sv_row = mm_sv.row(i);
          mm_sv_col = mm_sv.row(i);
          hess_sv -= mm_sv_col * dd * mm_sv_row;
        }
      }
    } else {
      hess_v = MatrixXd::Constant(form_len[0], form_len[0],
                                  std::numeric_limits<double>::quiet_NaN());
      hess_a = MatrixXd::Constant(form_len[1], form_len[1],
                                  std::numeric_limits<double>::quiet_NaN());
      hess_t0 = MatrixXd::Constant(form_len[2], form_len[2],
                                   std::numeric_limits<double>::quiet_NaN());
      hess_w = MatrixXd::Constant(form_len[3], form_len[3],
                                  std::numeric_limits<double>::quiet_NaN());
      hess_sv = MatrixXd::Constant(form_len[4], form_len[4],
                                   std::numeric_limits<double>::quiet_NaN());
    }
  }

  // return an R list (via Rcpp::List) of the Hessians
  return List::create(Named("drift") = -hess_v,
                      Named("boundary") = -hess_a,
                      Named("ndt") = -hess_t0,
                      Named("bias") = -hess_w,
                      Named("sv") = -hess_sv);
}


// Variance-Covariance Matrix Function
List fddm_fit::calc_vcov()
{
  // check if Hessians have been calculated
  if (hess_v.rows() < 1 && hess_a.rows() < 1 && hess_t0.rows() < 1 &&
      hess_w.rows() < 1 && hess_sv.rows() < 1) {
    calc_hessians(coefs);
  }

  // invert (already negated) hessians to get variance-covariance matrices
  vcov_v = hess_v.inverse();
  vcov_a = hess_a.inverse();
  vcov_t0 = hess_t0.inverse();
  vcov_w = hess_w.inverse();
  vcov_sv = hess_sv.inverse();

  // return an R list (via Rcpp::List) of the variance-covariance matrices
  return List::create(Named("drift") = vcov_v,
                      Named("boundary") = vcov_a,
                      Named("ndt") = vcov_t0,
                      Named("bias") = vcov_w,
                      Named("sv") = vcov_sv);
}


// Standard Error Function
VectorXd fddm_fit::calc_std_err()
{
  // check if variance-covariance matrices have been calculated
  if (vcov_v.rows() < 1 && vcov_a.rows() < 1 && vcov_t0.rows() < 1 &&
      vcov_w.rows() < 1 && vcov_sv.rows() < 1) {
    calc_vcov();
  }

  // calculate standard errors
  VectorXd std_err(Ncoefs);
  int par_idx = 0;
  if (form_len[0] > 0) {
    std_err.segment(par_idx, form_len[0]) = vcov_v.diagonal().cwiseSqrt();
    par_idx += form_len[0];
  }
  if (form_len[1] > 0) {
    std_err.segment(par_idx, form_len[1]) = vcov_a.diagonal().cwiseSqrt();
    par_idx += form_len[1];
  }
  if (form_len[2] > 0) {
    std_err.segment(par_idx, form_len[2]) = vcov_t0.diagonal().cwiseSqrt();
    par_idx += form_len[2];
  }
  if (form_len[3] > 0) {
    std_err.segment(par_idx, form_len[3]) = vcov_w.diagonal().cwiseSqrt();
    par_idx += form_len[3];
  }
  if (form_len[4] > 0) {
    std_err.segment(par_idx, form_len[4]) = vcov_sv.diagonal().cwiseSqrt();
    par_idx += form_len[4];
  }

  return std_err;
}
