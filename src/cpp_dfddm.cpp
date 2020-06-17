// Density function for the Ratcliff Diffusion Decision Model (DDM) PDF

#include "funcs.h"

using Rcpp::stop;
using std::vector;
using std::string;




// [[Rcpp::export]]
NumericVector cpp_dfddm(const NumericVector& rt,
                        const SEXP& response,
                        const NumericVector& a, const NumericVector& v,
                        const NumericVector& t0, const NumericVector& w,
                        const NumericVector& sv, const LogicalVector& log_prob,
                        const std::string& n_terms_small,
                        const std::string& summation_small,
                        const std::string& scale,
                        const NumericVector& eps)
{
  // convert responses to false (0, lower) and true (1, upper)
  // add factor
  vector<bool> resp;
  int Nres, type = TYPEOF(response);
  if (type == 10) { // LogicalVector
    resp = Rcpp::as<vector<bool> >(response);
    Nres = resp.size();
  } else if (type == 13 || type == 14) { // IntegerVector or NumericVector
    vector<double> temp = Rcpp::as<vector<double> >(response);
    Nres = temp.size();
    resp.reserve(Nres);
    double m = temp[0];
    for (int i = 1; i < Nres; i++) { // get min
      if (temp[i] < m) {
        m = temp[i];
      }
    }
    for (int i = 0; i < Nres; i++) {
      if (temp[i] > m) {
        resp[i] = 1;
      } else {
        resp[i] = 0;
      }
    }
  } else if (type == 16) { // StringVector (contains at least one string)
    vector<string> temp = Rcpp::as<vector<string> >(response);
    Nres = temp.size();
    resp.reserve(Nres);
    for (int i = 0; i < Nres; i++) {
      if (temp[i][0] == 'u' || temp[i][0] == 'U') {
        resp[i] = 1;
      } else if (temp[i][0] == 'l' || temp[i][0] == 'L') {
        resp[i] = 0;
      } else {
        stop("dfddm error: unknown value in function parameter 'response' at index %i.", i+1);
      }
    }
  } else {
    stop("dfddm error: type of function parameter 'response' vector is not one of: integer, double, boolean (logical), or string (character).");
  }

  // find Nmax (max length of parameter inputs)
  int Nrt  = rt.length();
  int Na   = a.length();
  int Nv   = v.length();
  int Nt0  = t0.length();
  int Nw   = w.length();
  int Nsv  = sv.length();
  int Nlog = log_prob.length();
  int Neps = eps.length();
  int Nmax = max({Nrt, Nres, Na, Nv, Nt0, Nw, Nsv, Nlog, Neps});

  // input checking
  for (int i = 0; i < Na; i++) {
    if (a[i] < 0) {
      stop("dfddm error: model parameter 'a' < 0 at index %i.", i+1);
    }
  }
  for (int i = 0; i < Nt0; i++) {
    if (t0[i] < 0) {
      stop("dfddm error: model parameter 't0' < 0 at index %i.", i+1);
    }
  }
  for (int i = 0; i < Nw; i++) {
    if (w[i] < 0 || w[i] > 1) {
      stop("dfddm error: model parameter 'w' < 0 or 'w' > 1 at index %i.", i+1);
    }
  }
  if (sv[0] != -1) {
    for (int i = 0; i < Nsv; i++) {
      if (sv[i] < 0) {
        stop("dfddm error: model parameter 'sv' < 0 at index %i.", i+1);
      }
    }
  }



  // prepare output
  NumericVector out(Nmax);
  double t;

  // determine which method to use
  char n_terms_small0 = n_terms_small[0];
  char summation_small0 = summation_small[summation_small.length()-1];
  SummFunc summ;

  if (n_terms_small0 == 'F' || n_terms_small0 == 'f') { // Foster method
    if (summation_small0 == '7') {
      summ = &small_sum_eps_17;
    } else if (summation_small0 == '4') {
      summ = &small_sum_eps_14;
    } else {
      stop("dfddm error: invalid function parameter summation_small");
      return NAN;
    }

    // loop through all inputs
    for (int i = 0; i < Nmax; i++) {
      t = rt[i % Nrt] - t0[i % Nt0]; // take non-decisison time from response time
      if (t <= 0) { // handle density outside of time bounds
        if (log_prob[i % Nlog]) {
          out[i] = -std::numeric_limits<double>::infinity();
        } else {
          out[i] = 0;
        }
        continue;
      }

      if (resp[i % Nres]) { // response is "upper" so use alternate parameters
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
        stop("dfddm error: invalid function parameter scale");
        return NAN;
      }
      if (n_terms_small0 == 'K' || n_terms_small0 == 'k') { // Kesselmeier
        numm = &ks_Kes;
      } else if (n_terms_small0 == 'N' || n_terms_small0 == 'n') { // Navarro
        numm = &ks_Nav;
      } else {
        stop("dfddm error: invalid function parameter n_terms_small");
        return NAN;
      }
      if (summation_small0 == '7') {
        summ = &small_sum_2017;
      } else if (summation_small0 == '4') {
        summ = &small_sum_2014;
      } else {
        stop("dfddm error: invalid function parameter summation_small");
        return NAN;
      }
    }

    // loop through all inputs
    for (int i = 0; i < Nmax; i++) {
      t = rt[i % Nrt] - t0[i % Nt0]; // take non-decisison time from response time
      if (t <= 0) { // handle density outside of time bounds
        if (log_prob[i % Nlog]) {
          out[i] = -std::numeric_limits<double>::infinity();
        } else {
          out[i] = 0;
        }
        continue;
      }

      if (resp[i % Nres]) { // response is "upper" so use alternate parameters
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
