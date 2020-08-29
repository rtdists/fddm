// Density function for the Ratcliff Diffusion Decision Model (DDM) PDF

#include "funcs.h"

using Rcpp::warning;
using Rcpp::stop;
using std::vector;
using std::string;




// [[Rcpp::export]]
NumericVector cpp_dfddm(const NumericVector& rt,
                        const SEXP& response,
                        const NumericVector& a, const NumericVector& v,
                        const NumericVector& t0, const NumericVector& w,
                        const NumericVector& sv, const bool& log_prob,
                        const std::string& n_terms_small,
                        const std::string& summation_small,
                        const std::string& scale,
                        const int& max_terms_large,
                        const NumericVector& eps)
{
  // convert responses to false (0, lower) and true (1, upper)
  vector<bool> resp;
  int Nres, type = TYPEOF(response);
  if (type == 10) { // LogicalVector
    resp = Rcpp::as<vector<bool> >(response);
    Nres = resp.size();
  } else if (type == 13) { // IntegerVector (including factors)
    Rcpp::IntegerVector temp = Rcpp::wrap(response);
    if (Rf_isFactor(response) == 1) { // factor
      Rcpp::CharacterVector levs = temp.attr("levels");
      vector<int> temp = Rcpp::as<vector<int> >(response);
      Nres = temp.size();
      resp.reserve(Nres);
      for (int i = 0; i < Nres; i++) {
        if (levs[temp[i]-1] == levs[0]) { // lower
          resp[i] = 0;
        } else if (levs[temp[i]-1] == levs[1]) { // upper
          resp[i] = 1;
        } else {
          stop("dfddm error: index %i of function parameter 'response' contains a value that is neither 1 nor 2", i+1);
        }
      }
    } else { // IntegerVector, NOT factor
      vector<int> temp = Rcpp::as<vector<int> >(response);
      Nres = temp.size();
      resp.reserve(Nres);
      for (int i = 0; i < Nres; i++) {
        if (temp[i] == 1) { // lower
          resp[i] = 0;
        } else if (temp[i] == 2){ // upper
          resp[i] = 1;
        } else {
          stop("dfdmm error: function parameter 'response' was input as a vector of integers, and an integer other than 1 or 2 was detected at index %i", i+1);
        }
      }
    }
  } else if (type == 14) { // NumericVector
    vector<double> temp = Rcpp::as<vector<double> >(response);
    Nres = temp.size();
    resp.reserve(Nres);
    for (int i = 0; i < Nres; i++) {
      if (temp[i] == 1) { // lower
        resp[i] = 0;
      } else if (temp[i] == 2){ // upper
        resp[i] = 1;
      } else {
        stop("dfdmm error: function parameter 'response' was input as a vector of integers, and an integer other than 1 or 2 was detected at index %i", i+1);
      }
    }
  } else if (type == 16) { // StringVector (contains at least one string)
    vector<string> temp = Rcpp::as<vector<string> >(response);
    Nres = temp.size();
    resp.reserve(Nres);
    for (int i = 0; i < Nres; i++) {
      if (temp[i][0] == 'l' || temp[i][0] == 'L') { // lower
        resp[i] = 0;
      } else if (temp[i][0] == 'u' || temp[i][0] == 'U') { // upper
        resp[i] = 1;
      } else {
        stop("dfddm error: function parameter 'response' was input as a vector of strings (characters), and an object other than 'u' or 'l' (case insensitive) was detected as the first character at index %i.", i+1);
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
  int Neps = eps.length();
  int Nmax = max({Nrt, Nres, Na, Nv, Nt0, Nw, Nsv, Neps});

  // input checking
  if (Nrt < 1) {
    warning("dfddm warning: model input 'rt' is empty, empty vector returned\n");
    NumericVector empty_out(0);
    return empty_out;
  }
  if (Nres < 1) {
    stop("dfddm error: model input 'response' is empty");
  }
  if (Na < 1) {
    stop("dfddm error: model parameter 'a' is empty");
  }
  if (Nv < 1) {
    stop("dfddm error: model parameter 'v' is empty");
  }
  if (Nt0 < 1) {
    stop("dfddm error: model parameter 't0' is empty");
  }
  if (Nw < 1) {
    stop("dfddm error: model parameter 'w' is empty");
  }
  if (Nsv < 1) {
    stop("dfddm error: model parameter 'sv' is empty");
  }
  if (Neps < 1) {
    stop("dfddm error: model parameter 'err_tol' is empty");
  }
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



  // determine which method to use
  char n_terms_small0 = n_terms_small[0];
  char summation_small0 = summation_small[summation_small.length()-1];
  char scale0 = scale[0];
  NummFunc numm;
  SummFunc summ;
  DensFunc dens;
  double rt0;

  if (log_prob) { // calculate log(probability)
    rt0 = -std::numeric_limits<double>::infinity();
    if (n_terms_small0 == 'S' || n_terms_small0 == 's') { // SWSE method
      if (scale0 == 'b' || scale0 == 'B') { // both
        dens = &fc_log;
      } else if (scale0 == 's' || scale0 == 'S'){ // small
        dens = &ff_log;
      } else {
        stop("dfddm error: invalid function parameter 'scale': %s", scale);
      }
      numm = NULL;
      if (summation_small0 == '7') { // 2017
        summ = &small_sum_eps_17;
      } else if (summation_small0 == '4') { // 2014
        summ = &small_sum_eps_14;
      } else {
        stop("dfddm error: invalid function parameter 'summation_small': %s",
             summation_small);
      }
    } else {
      if (scale0 == 'l' || scale0 == 'L') { // large
        dens = &fl_log;
        numm = NULL;
        summ = NULL;
      } else {
        if (scale0 == 'b' || scale0 == 'B') { // both
          dens = &fb_log;
        } else if (scale0 == 's' || scale0 == 'S') { // small
          dens = &fs_log;
        } else {
          stop("dfddm error: invalid function parameter 'scale': %s", scale);
        }
        if (n_terms_small0 == 'G' || n_terms_small0 == 'g') { // Gondan
          numm = &ks_Gon;
        } else if (n_terms_small0 == 'N' || n_terms_small0 == 'n') { // Navarro
          numm = &ks_Nav;
        } else {
          stop("dfddm error: invalid function parameter 'n_terms_small': %s",
               n_terms_small);
        }
        if (summation_small0 == '7') { // 2017
          summ = &small_sum_2017;
        } else if (summation_small0 == '4') { // 2014
          summ = &small_sum_2014;
        } else {
          stop("dfddm error: invalid function parameter 'summation_small': %s",
               summation_small);
        }
      }
    }
  } else { // calculate regular probability
    rt0 = 0;
    if (n_terms_small0 == 'S' || n_terms_small0 == 's') { // SWSE method
      if (scale0 == 'b' || scale0 == 'B') { // both
        dens = &fc;
      } else if (scale0 == 's' || scale0 == 'S'){ // small
        dens = &ff;
      } else {
        stop("dfddm error: invalid function parameter 'scale': %s", scale);
      }
      numm = NULL;
      if (summation_small0 == '7') { // 2017
        summ = &small_sum_eps_17;
      } else if (summation_small0 == '4') { // 2014
        summ = &small_sum_eps_14;
      } else {
        stop("dfddm error: invalid function parameter 'summation_small': %s",
             summation_small);
      }
    } else {
      if (scale0 == 'l' || scale0 == 'L') { // large
        dens = &fl;
        numm = NULL;
        summ = NULL;
      } else {
        if (scale0 == 'b' || scale0 == 'B') { // both
          dens = &fb;
        } else if (scale0 == 's' || scale0 == 'S') { // small
          dens = &fs;
        } else {
          stop("dfddm error: invalid function parameter 'scale': %s", scale);
        }
        if (n_terms_small0 == 'G' || n_terms_small0 == 'g') { // Gondan
          numm = &ks_Gon;
        } else if (n_terms_small0 == 'N' || n_terms_small0 == 'n') { // Navarro
          numm = &ks_Nav;
        } else {
          stop("dfddm error: invalid function parameter 'n_terms_small': %s",
               n_terms_small);
        }
        if (summation_small0 == '7') { // 2017
          summ = &small_sum_2017;
        } else if (summation_small0 == '4') { // 2014
          summ = &small_sum_2014;
        } else {
          stop("dfddm error: invalid function parameter 'summation_small': %s",
               summation_small);
        }
      }
    }
  }



  // loop through all inputs
  NumericVector out(Nmax);
  double t;
  for (int i = 0; i < Nmax; i++) {
    t = rt[i % Nrt] - t0[i % Nt0]; // take non-decisison time from response time
    if (t <= 0) { // handle density outside of time bounds
      out[i] = rt0;
      continue;
    }
    if (resp[i % Nres]) { // response is "upper" so use alternate parameters
      out[i] = dens(t, a[i % Na], -v[i % Nv], 1 - w[i % Nw], sv[i % Nsv],
                    eps[i % Neps], max_terms_large, numm, summ);
    } else { // response is "lower" so use unchanged parameters
      out[i] = dens(t, a[i % Na], v[i % Nv], w[i % Nw], sv[i % Nsv],
                    eps[i % Neps], max_terms_large, numm, summ);
    }
  }


  return out;
}
