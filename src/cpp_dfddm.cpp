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
  
  if (scale == "small") {
    if (summation_small == "2017") {
      if (n_terms_small == "Foster") {
        dens = fs_Fos_17;
      } else if (n_terms_small == "Kesselmeier") {
        dens = fs_Kes_17;
      } else if (n_terms_small == "Navarro") {
        dens = fs_Nav_17;
      } else {
        Rcerr << "error: invalid n_terms_small" << endl;
        return NAN;
      }
    } else if (summation_small == "2014") {
      if (n_terms_small == "Foster") {
        dens = fs_Fos_14;
      } else if (n_terms_small == "Kesselmeier") {
        dens = fs_Kes_14;
      } else if (n_terms_small == "Navarro") {
        dens = fs_Nav_14;
      } else {
        Rcerr << "error: invalid n_terms_small" << endl;
        return NAN;
      }
    } else {
      Rcerr << "error: invalid summation_small" << endl;
      return NAN;
    }
  } else if (scale == "large") {
    dens = fl_Nav_09;
  } else if (scale == "both") {
    if (summation_small == "2017") {
      if (n_terms_small == "Kesselmeier") {
        dens = fb_Kes_17;
      } else if (n_terms_small == "Navarro") {
        dens = fb_Nav_17;
      } else {
        Rcerr << "error: invalid n_terms_small" << endl;
        return NAN;
      }
    } else if (summation_small == "2014") {
      if (n_terms_small == "Kesselmeier") {
        dens = fb_Kes_14;
      } else if (n_terms_small == "Navarro") {
        dens = fb_Nav_14;
      } else {
        Rcerr << "error: invalid n_terms_small" << endl;
        return NAN;
      }
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
      out[i] = dens(t, a[i % Na], -v[i % Nv], 1 - w[i % Nw],
                    sv[i % Nsv], log_prob[i % Nlog], eps[i % Neps]);
    } else { // response is "lower" so use unchanged parameters
      out[i] = dens(t, a[i % Na], v[i % Nv], w[i % Nw],
                    sv[i % Nsv], log_prob[i % Nlog], eps[i % Neps]);
    }
    
  }

  return out;
}


