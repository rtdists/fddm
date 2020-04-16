// Density function for the Ratcliff Diffusion Decision Model (DDM) PDF

#include "funcs.hpp"

using std::endl;
using Rcpp::Rcerr;




// [[Rcpp::export]]
NumericVector cpp_dfddm_fast(const NumericVector& rt,
                             const double& a, const double& v,
                             const double& t0, const double& w,
                             const double& eps)
{
  int Nmax = rt.length();
  NumericVector out(Nmax);
  double t;
  double mult_s = a * exp(-v * a * w) / SQRT_2PI;

  if (eps >= 1e-6) {
  double mult_l = M_PI * exp(-v * a * w) / (a*a);
  double gamma = -M_PI*M_PI / (2 * a*a);
    for (int i = 0; i < Nmax; i++) {
      t = rt[i] - t0;
      if (t >= 2.5) { // use large-time (1 term in summation)
        out[i] = mult_l * exp(-v*v * t / 2) * sin(w * M_PI) * exp(gamma*t);
      } else { // use small-time
        out[i] = mult_s * exp(-v * v * t / 2) / (t * sqrt(t)) *
                  small_sum_eps_17(t, a, w, 0, eps / mult_s);
      }

    }
  } else {
    Rcerr << eps << " < 1e-6" << endl;
    for (int i = 0; i < Nmax; i++) {
      t = rt[i] - t0;
      out[i] = mult_s * exp(-v * v * t / 2) / (t * sqrt(t)) *
                small_sum_eps_17(t, a, w, 0, eps / mult_s);
    }
  }
  return out;
}


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

  // determine which method to use (and parameter validation)
  char n_terms_small0 = n_terms_small[0];
  char summation_small0 = summation_small[summation_small.length()-1];
  char scale0 = scale[0];

  DensFunc dens;
  SummFunc summ;
  NummFunc numm;

  if (n_terms_small0 == 'F' || n_terms_small0 == 'f') { // Foster
    dens = &fs; // density must be small-time
    numm = NULL;
    if (summation_small0 == '7') {
      summ = &small_sum_eps_17;
    } else if (summation_small0 == '4') {
      summ = &small_sum_eps_14;
    } else {
      Rcerr << "error: invalid parameter summation_small" << endl;
      return NAN;
    }
  } else {
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
  }

  // loop through all inputs
  NumericVector out(Nmax);
  double t;

  for (int i = 0; i < Nmax; i++) {
    t = rt[i % Nrt] - t0[i % Nt0]; // take non-decisison time from response time
    if (t <= 0) {
      out[i] = 0; // if t=0, return 0 instead of potential NAN
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
