// Density function for the Ratcliff Diffusion Decision Model (DDM) PDF

#include "funcs.h"

using std::endl;
using Rcpp::Rcerr;




// [[Rcpp::export]]
NumericVector cpp_dfddm(const NumericVector& rt,
                        const LogicalVector& response,
                        const NumericVector& a, const NumericVector& v,
                        const NumericVector& t0, const NumericVector& w,
                        const NumericVector& sv, const LogicalVector& log_prob,
                        const std::string& n_terms_small,
                        const std::string& summation_small,
                        const std::string& scale,
                        const NumericVector& eps)
{
  // find Nmax (max length of parameter inputs)
  int Nrt  = rt.length();
  int Nres = response.length();
  int Na   = a.length();
  int Nv   = v.length();
  int Nt0  = t0.length();
  int Nw   = w.length();
  int Nsv  = sv.length();
  int Nlog = log_prob.length();
  int Neps = eps.length();
  int Nmax = max({Nrt, Nres, Na, Nv, Nt0, Nw, Nsv, Nlog, Neps});
  
  // determine which method to use
  DensFunc dens;
  SummFunc summ;
  NummFunc numm;

  if (scale == "small") {
    dens = &fs;
    if (n_terms_small == "Foster") {
      numm = NULL;
      if (summation_small == "2017") {
        summ = &small_sum_eps_17;
      } else if (summation_small == "2014") {
        summ = &small_sum_eps_14;
      } else {
        Rcerr << "error: invalid summation_small" << endl;
        return NAN;
      }
    } else {
      if (n_terms_small == "Kesselmeier") {
        numm = &ks_Kes;
      } else if (n_terms_small == "Navarro") {
        numm = &ks_Nav;
      } else {
        Rcerr << "error: invalid n_terms_small" << endl;
        return NAN;
      }
      if (summation_small == "2017") {
        summ = &small_sum_2017;
      } else if (summation_small == "2014") {
        summ = &small_sum_2014;
      } else {
        Rcerr << "error: invalid summation_small" << endl;
        return NAN;
      }
    }
  } else if (scale == "large") {
    dens = &fl;
    numm = NULL;
    summ = NULL;
  } else if (scale == "both") {
    dens = &fb;
    if (n_terms_small == "Kesselmeier") {
      numm = &ks_Kes;
    } else if (n_terms_small == "Navarro") {
      numm = &ks_Nav;
    } else {
      Rcerr << "error: invalid n_terms_small" << endl;
      return NAN;
    }
    if (summation_small == "2017") {
      summ = &small_sum_2017;
    } else if (summation_small == "2014") {
      summ = &small_sum_2014;
    } else {
      Rcerr << "error: invalid summation_small" << endl;
      return NAN;
    }
  } else {
    Rcerr << "error: invalid scale" << endl;
    return NAN;
  }

  
  NumericVector out(Nmax);
  double t;
  
  for (int i = 0; i < Nmax; i++) {
    t = rt[i % Nrt] - t0[i % Nt0]; // take non-decisison time from response time
    if (t <= 0) {
      out[i] = 0; // if t=0, return 0 instead of NaN
      continue;
    }
    
    if (response[i % Nres]) { // response is "upper" so use alternate parameters
      out[i] = dens(t, a[i % Na], -v[i % Nv], 1 - w[i % Nw], sv[i % Nsv],
                    log_prob[i % Nlog], eps[i % Neps], numm, summ, -1);
    } else { // response is "lower" so use unchanged parameters
      out[i] = dens(t, a[i % Na], v[i % Nv], w[i % Nw], sv[i % Nsv],
                    log_prob[i % Nlog], eps[i % Neps], numm, summ, -1);
    }
    
  }

  return out;
}
