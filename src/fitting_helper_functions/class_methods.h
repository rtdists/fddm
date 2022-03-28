// included from src/fitting_class.cpp

#include "declarations.h" // all function declarations and constant definitions

#include "par_checks.h" // src/fitting_helper_functions/par_checks.h
#include "likelihood_funcs.h" // src/fitting_helper_functions/likelihood_funcs.h
#include "gradient_funcs.h" // src/fitting_helper_functions/gradient_funcs.h
#include "sum_funcs.h" // src/fitting_helper_functions/sum_funcs.h
#include "num_funcs.h" // src/fitting_helper_functions/num_funcs.h

using std::log;
using std::isfinite;



// Constructors
fddmFit::fddmFit()
{
  Rcpp::warning("fddmFit warning: empty constructor used, need to manually input the response times and responses.");
}

fddmFit::fddmFit(const vector<double>& rt_vector,
                   const SEXP& response_vector)
{
  rt = check_rt(rt_vector, Nrt); // also sets Nrt
  response = convert_responses(response_vector, Nrt);
  likelihood.resize(Nrt);
}

fddmFit::fddmFit(const vector<double>& rt_vector,
                   const SEXP& response_vector,
                   const double& error_tolerance)
{
  rt = check_rt(rt_vector, Nrt); // also sets Nrt
  response = convert_responses(response_vector, Nrt);
  likelihood.resize(Nrt);
  err_tol = check_err_tol(error_tolerance);
}

fddmFit::fddmFit(const vector<double>& rt_vector,
                   const SEXP& response_vector,
                   const double& error_tolerance,
                   const double& switching_threshold)
  {
  rt = check_rt(rt_vector, Nrt); // also sets Nrt
  response = convert_responses(response_vector, Nrt);
  likelihood.resize(Nrt);
  err_tol = check_err_tol(error_tolerance);
  switch_thresh = check_switch_thresh(switching_threshold);
}


// Log-likelihood Function
double fddmFit::calc_loglik(const vector<double>& temp_params)
{
  // check parameters and store them in the order: a, v, t0, w, sv
  if (invalid_parameters(parameters, temp_params)) {
    for (int i = 0; i < Nrt; i++) {
      likelihood[i] = rt0;
    }
    return rt0;
  }

  // calculate loglikelihood and -sum over all data
  double ll = 0.0;
  double t;
  for (int i = 0; i < Nrt; i++) {
    t = rt[i] - parameters[2];
    if (t > 0 && isfinite(t)) {
      if (response[i] == 1) { // response is "lower"
        likelihood[i] = ft(t, parameters[0], parameters[1], parameters[3],
                           parameters[4], err_tol, switch_thresh);
      } else { // response is "upper" so use alternate parameters
        likelihood[i] = ft(t, parameters[0], -parameters[1], 1-parameters[3],
                           parameters[4], err_tol, switch_thresh);
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
vector<double> fddmFit::calc_gradient(const vector<double>& temp_params)
{
  // check if new parameters match the stored ones
  if (temp_params != parameters) { // if not, must recalculate likelihood
    double not_used = calc_loglik(temp_params);
  } // if they match, can reuse the likelihood
  
  // calculate gradient and -sum over all data
  // note that d/dx(log(f(x))) = d/dx(f(x)) / f(x) (why dividing by likelihood)
  vector<double> gradient(5, 0.0);
  double t;
  for (int i = 0; i < Nrt; i++) {
    t = rt[i] - parameters[2];
    if (t > 0 && isfinite(t)) { // maybe don't have to check b/c same in lik calc?
      if (response[i] == 1) { // response is "lower"
        gradient[0] -= da(t, parameters[0], parameters[1], parameters[3],
                          parameters[4], err_tol, switch_thresh)
                       / likelihood[i];
        gradient[1] -= dv(t, parameters[0], parameters[1], parameters[3],
                          parameters[4], err_tol, switch_thresh)
                       / likelihood[i];
        gradient[2] -= dt0(t, parameters[0], parameters[1], parameters[3],
                           parameters[4], err_tol, switch_thresh)
                       / likelihood[i];
        gradient[3] -= dw(t, parameters[0], parameters[1], parameters[3],
                          parameters[4], err_tol, switch_thresh)
                       / likelihood[i];
        gradient[4] -= dsv(t, parameters[0], parameters[1], parameters[3],
                           parameters[4], err_tol, switch_thresh)
                       / likelihood[i];
      } else { // response is "upper" so use alternate parameters
        gradient[0] -= da(t, parameters[0], -parameters[1], 1-parameters[3],
                          parameters[4], err_tol, switch_thresh)
                       / likelihood[i];
        gradient[1] += dv(t, parameters[0], -parameters[1], 1-parameters[3],
                          parameters[4], err_tol, switch_thresh)
                       / likelihood[i]; // chain rule negates derivative
        gradient[2] -= dt0(t, parameters[0], -parameters[1], 1-parameters[3],
                           parameters[4], err_tol, switch_thresh)
                       / likelihood[i];
        gradient[3] += dw(t, parameters[0], -parameters[1], 1-parameters[3],
                          parameters[4], err_tol, switch_thresh)
                       / likelihood[i]; // chain rule negates derivative
        gradient[4] -= dsv(t, parameters[0], -parameters[1], 1-parameters[3],
                           parameters[4], err_tol, switch_thresh)
                       / likelihood[i];
      }
    } else {
      likelihood[i] = rt0;
      for (int j = 0; j < 5; j++) {
        gradient[j] = std::numeric_limits<double>::quiet_NaN();
      }
      break;
    }
  }
  return gradient;
}
