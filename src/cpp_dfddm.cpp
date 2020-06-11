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

  // prepare output
  NumericVector out(Nmax);
  double t;

  // determine which method to use (and parameter validation)
  // paramater validation idea:
    // if the input is invalid, should we default to something?
  char n_terms_small0 = n_terms_small[0];
  char summation_small0 = summation_small[summation_small.length()-1];
  SummFunc summ;

  if (n_terms_small0 == 'F' || n_terms_small0 == 'f') { // Foster method
    if (summation_small0 == '7') {
      summ = &small_sum_eps_17;
    } else if (summation_small0 == '4') {
      summ = &small_sum_eps_14;
    } else {
      Rcerr << "error: invalid parameter summation_small" << endl;
      return NAN;
    }

    // loop through all inputs
    for (int i = 0; i < Nmax; i++) {
      t = rt[i % Nrt] - t0[i % Nt0]; // take non-decisison time from response time
      if (t <= 0) {
        out[i] = 0; // if t=0, return 0 instead of potential NAN
        continue;
      }

      if (response[i % Nres]) { // response is "upper" so use alternate parameters
        out[i] = ff(t, a[i % Na], -v[i % Nv], 1 - w[i % Nw], sv[i % Nsv],
                      log_prob[i % Nlog], eps[i % Neps], NULL, summ);
      } else { // response is "lower" so use unchanged parameters
        out[i] = ff(t, a[i % Na], v[i % Nv], w[i % Nw], sv[i % Nsv],
                      log_prob[i % Nlog], eps[i % Neps], NULL, summ);
      }
    }



  } else{ // Kesselmeier and Navarro methods
    char scale0 = scale[0];
    DensFunc dens;
    NummFunc numm;
    if (scale0 == 'l' || scale0 == 'L') {
      dens = &fl;
      numm = NULL;
      summ = NULL;
    } else {
      if (scale0 == 's' || scale0 == 'S') {
        dens = &fs;
      } else if (scale0 == 'b' || scale0 == 'B') {
        dens = &fb;
      } else {
        Rcerr << "error: invalid parameter scale" << endl;
        return NAN;
      }
      if (n_terms_small0 == 'K' || n_terms_small0 == 'k') { // Kesselmeier
        numm = &ks_Kes;
      } else if (n_terms_small0 == 'N' || n_terms_small0 == 'n') { // Navarro
        numm = &ks_Nav;
      } else {
        Rcerr << "error: invalid parameter n_terms_small" << endl;
        return NAN;
      }
      if (summation_small0 == '7') {
        summ = &small_sum_2017;
      } else if (summation_small0 == '4') {
        summ = &small_sum_2014;
      } else {
        Rcerr << "error: invalid parameter summation_small" << endl;
        return NAN;
      }
    }

    // loop through all inputs
    for (int i = 0; i < Nmax; i++) {
      t = rt[i % Nrt] - t0[i % Nt0]; // take non-decisison time from response time
      if (t <= 0) {
        out[i] = 0; // if t=0, return 0 instead of potential NAN
        continue;
      }

      if (response[i % Nres]) { // response is "upper" so use alternate parameters
        out[i] = dens(t, a[i % Na], -v[i % Nv], 1 - w[i % Nw], sv[i % Nsv],
                      log_prob[i % Nlog], eps[i % Neps], numm, summ);
      } else { // response is "lower" so use unchanged parameters
        out[i] = dens(t, a[i % Na], v[i % Nv], w[i % Nw], sv[i % Nsv],
                      log_prob[i % Nlog], eps[i % Neps], numm, summ);
      }
    }
  }



  return out;
}
