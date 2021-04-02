// Helper functions for cpp_dfddm function (located at src/cpp_dfddm.cpp)

#include "funcs.h"

vector<int> convert_responses(const SEXP& response, int& Nres)
{
  vector<int> resp, bad_idx;
  int type = TYPEOF(response);
  int bad_par = 0;

  if (type == 13 || type == 14) { // IntegerVector (including factors) or NumericVector
    resp = Rcpp::as<vector<int> >(response);
    Nres = resp.size();
    for (int i = 0; i < Nres; i++) {
      if (resp[i] != 1 && resp[i] != 2) { // factors automatically sort themselves
        resp[i] = 0;
        bad_par++;
        bad_idx.push_back(i);
      }
    }
  } else if (type == 16) { // StringVector (contains at least one string)
    vector<string> temp = Rcpp::as<vector<string> >(response);
    Nres = temp.size();
    for (int i = 0; i < Nres; i++) {
      if (temp[i][0] == 'l' || temp[i][0] == 'L') { // lower
        resp.push_back(1);
      } else if (temp[i][0] == 'u' || temp[i][0] == 'U') { // upper
        resp.push_back(2);
      } else {
        resp.push_back(0);
        bad_par++;
        bad_idx.push_back(i);
      }
    }
  } else if (type == 10) { // LogicalVector (boolean values)
    resp = Rcpp::as<vector<int> >(response);
    Nres = resp.size();
    for (int i = 0; i < Nres; i++) {
      if (resp[i] == 0) { // lower
        resp[i] = 1;
      } else if (resp[i] == 1){ // upper
        resp[i] = 2;
      } else {
        resp[i] = 0;
        bad_par++;
        bad_idx.push_back(i);
      }
    }
  } else {
    stop("dfddm error: type of function parameter 'response' vector is not one of: integer, double, factor, string (character), or boolean (logical).");
  }

  if (bad_par > 0) { // error handling
    std::string warning_text = "dfddm warning: function parameter 'response' was input as a vector of ";

    if (type == 13) { // IntegerVector
      if (Rf_isFactor(response) == 1) { // factor
        warning_text = warning_text.append("factors, and a value that does not match either the first or second level");
      } else { // NOT factor
        warning_text = warning_text.append("integers, and a value other than 1 or 2");
      }
    } else if (type == 14) {
      warning_text = warning_text.append("doubles (truncated to integers), and a value other than 1 or 2");
    } else if (type == 16) {
      warning_text = warning_text.append("strings (characters), and a value other than 'L' or 'U' (case insensitive)");
    } else if (type == 10) {
      warning_text = warning_text.append("booleans (logicals), and a value other than TRUE or FALSE was");
    } else {
      warning_text = warning_text.append("some unknown type");
    }

    warning_text = warning_text.append(" was detected at ");

    if (bad_par == 1) {
      warning_text = warning_text.append("index ");
      warning_text = warning_text.append(to_string(bad_idx[0]+1));
      warning_text = warning_text.append("; zeros produced.");
      warning(warning_text);
    } else {
      warning_text = warning_text.append("the following indices: ");
      warning_text = warning_text.append(to_string(bad_idx[0]+1));
      for(int j = 1; j < bad_par; j++) {
        warning_text = warning_text.append(", ");
        warning_text = warning_text.append(to_string(bad_idx[j]+1));
      }
      warning_text = warning_text.append("; zeros produced.");
      warning(warning_text);
    }
  }

  return resp;
}



bool parameter_check(const int& Nrt, const int& Nres, const int& Na,
                     const int& Nv, const int& Nt0, const int& Nw,
                     const int& Nsv, const int& Nsig, const int& Nerr,
                     const int& Nmax, const NumericVector& rt,
                     const NumericVector& a, const NumericVector& v,
                     const NumericVector& t0, const NumericVector& w,
                     const NumericVector& sv, const NumericVector& sigma,
                     const NumericVector& err, vector<bool>& invalid_input)
{
  bool out = 1;
  if (Nrt < 1) {
    warning("dfddm warning: function parameter 'rt' is empty; empty vector returned.");
    out = 0;
  } // invalid inputs handled in calculate_pdf()
  if (Nres < 1) {
    warning("dfddm warning: function parameter 'response' is empty; empty vector returned.");
    out = 0;
  } // invalid inputs handled in calculate_pdf()
  if (Na < 1) {
    warning("dfddm warning: model parameter 'a' is empty; empty vector returned.");
    out = 0;
  } else {
    for (int i = 0; i < Na; i++) {
      if (a[i] <= 0) {
        for (int j = i; j < Nmax; j += Na) {
          invalid_input[j] = 1;
        }
      }
    }
  }
  if (Nv < 1) {
    warning("dfddm warning: model parameter 'v' is empty; empty vector returned.");
    out = 0;
  } else {
    for (int i = 0; i < Nv; i++) {
      if (std::isnan(v[i])) {
        for (int j = i; j < Nmax; j += Nv) {
          invalid_input[j] = 1;
        }
      }
    }
  }
  if (Nt0 < 1) {
    out = 0;
  } else {
    for (int i = 0; i < Nt0; i++) {
      if (t0[i] < 0) {
        for (int j = i; j < Nmax; j += Nt0) {
          invalid_input[j] = 1;
        }
      }
    }
  }
  if (Nw < 1) {
    warning("dfddm warning: model parameter 'w' is empty; empty vector returned.");
    out = 0;
  } else {
    for (int i = 0; i < Nw; i++) {
      if (w[i] <= 0 || w[i] >= 1) {
        for (int j = i; j < Nmax; j += Nw) {
          invalid_input[j] = 1;
        }
      }
    }
  }
  if (Nsv < 1) {
    warning("dfddm warning: model parameter 'sv' is empty; empty vector returned.");
    out = 0;
  } else {
    for (int i = 0; i < Nsv; i++) {
      if (sv[i] < 0) {
        for (int j = i; j < Nmax; j += Nsv) {
          invalid_input[j] = 1;
        }
      }
    }
  }
  if (Nsig < 1) {
    warning("dfddm warning: model parameter 'sigma' is empty; empty vector returned.");
    out = 0;
  } else {
    for (int i = 0; i < Nsig; i++) {
      if (sigma[i] <= 0) {
        for (int j = i; j < Nmax; j += Nsig) {
          invalid_input[j] = 1;
        }
      }
    }
  }
  if (Nerr < 1) {
    stop("dfddm error: function parameter 'err_tol' is empty.");
    out = 0;
  } else {
    bool bad_par = 0;
    vector<int> bad_idx;
    for (int i = 0; i < Nerr; i++) {
      if (err[i] <= 0) {
        bad_par = 1;
        bad_idx.push_back(i);
      }
    }
    if (bad_par) {
      int n_bad = bad_idx.size();
      if (n_bad == 1) {
        stop("dfddm error: function parameter 'err_tol' <= 0 at index %i.",
             bad_idx[0]+1);
      } else {
        std::string error_text = "dfddm error: function parameter 'err_tol' <= 0 at the following indices:\n";
        error_text = error_text.append(to_string(bad_idx[0]+1));
        for(int j = 1; j < n_bad; j++) {
          error_text = error_text.append(", ");
          error_text = error_text.append(to_string(bad_idx[j]+1));
        }
        error_text = error_text.append(".");
        stop(error_text);
      }
    }
  }

  return out;
}




void determine_method(const std::string& n_terms_small,
                      const std::string& summation_small,
                      const std::string& scale,
                      NumFunc& numf, SumFunc& sumf, DenFunc& denf,
                      double& rt0, const bool& log_prob)
{
  char n_terms_small0 = (!n_terms_small.empty()) ? n_terms_small[0] : EMPTYCHAR;
  char summation_small0 = (!summation_small.empty()) ?
    summation_small[summation_small.length()-1] : EMPTYCHAR;
  char scale0 = (!scale.empty()) ? scale[0] : EMPTYCHAR;

  if (log_prob) { // calculate log(probability)
    rt0 = -std::numeric_limits<double>::infinity();
    if (n_terms_small0 == 'S' || n_terms_small0 == 's') { // SWSE method
      if (scale0 == 'b' || scale0 == 'B') { // both
        denf = &fc_log;
      } else if (scale0 == 's' || scale0 == 'S'){ // small
        denf = &ff_log;
      } else {
        stop("dfddm error: invalid function parameter 'scale': %s.", scale);
      }
      numf = NULL;
      if (summation_small0 == '7') { // 2017
        sumf = &small_sum_eps_17;
      } else if (summation_small0 == '4') { // 2014
        sumf = &small_sum_eps_14;
      } else {
        stop("dfddm error: invalid function parameter 'summation_small': %s.",
             summation_small);
      }
    } else {
      if (scale0 == 'l' || scale0 == 'L') { // large
        denf = &fl_log;
        numf = NULL;
        sumf = NULL;
      } else {
        if (scale0 == 'b' || scale0 == 'B') { // both
          denf = &fb_log;
        } else if (scale0 == 's' || scale0 == 'S') { // small
          denf = &fs_log;
        } else {
          stop("dfddm error: invalid function parameter 'scale': %s.", scale);
        }
        if (n_terms_small0 == 'G' || n_terms_small0 == 'g') { // Gondan
          numf = &ks_Gon;
        } else if (n_terms_small0 == 'N' || n_terms_small0 == 'n') { // Navarro
          numf = &ks_Nav;
        } else {
          stop("dfddm error: invalid function parameter 'n_terms_small': %s.",
               n_terms_small);
        }
        if (summation_small0 == '7') { // 2017
          sumf = &small_sum_2017;
        } else if (summation_small0 == '4') { // 2014
          sumf = &small_sum_2014;
        } else {
          stop("dfddm error: invalid function parameter 'summation_small': %s.",
               summation_small);
        }
      }
    }
  } else { // calculate regular probability
    rt0 = 0;
    if (n_terms_small0 == 'S' || n_terms_small0 == 's') { // SWSE method
      if (scale0 == 'b' || scale0 == 'B') { // both
        denf = &fc;
      } else if (scale0 == 's' || scale0 == 'S'){ // small
        denf = &ff;
      } else {
        stop("dfddm error: invalid function parameter 'scale': %s.", scale);
      }
      numf = NULL;
      if (summation_small0 == '7') { // 2017
        sumf = &small_sum_eps_17;
      } else if (summation_small0 == '4') { // 2014
        sumf = &small_sum_eps_14;
      } else {
        stop("dfddm error: invalid function parameter 'summation_small': %s.",
             summation_small);
      }
    } else {
      if (scale0 == 'l' || scale0 == 'L') { // large
        denf = &fl;
        numf = NULL;
        sumf = NULL;
      } else {
        if (scale0 == 'b' || scale0 == 'B') { // both
          denf = &fb;
        } else if (scale0 == 's' || scale0 == 'S') { // small
          denf = &fs;
        } else {
          stop("dfddm error: invalid function parameter 'scale': %s.", scale);
        }
        if (n_terms_small0 == 'G' || n_terms_small0 == 'g') { // Gondan
          numf = &ks_Gon;
        } else if (n_terms_small0 == 'N' || n_terms_small0 == 'n') { // Navarro
          numf = &ks_Nav;
        } else {
          numf = &ks_Gon;
          stop("dfddm error: invalid function parameter 'n_terms_small': %s.",
               n_terms_small);
        }
        if (summation_small0 == '7') { // 2017
          sumf = &small_sum_2017;
        } else if (summation_small0 == '4') { // 2014
          sumf = &small_sum_2014;
        } else {
          sumf = &small_sum_2017;
          stop("dfddm error: invalid function parameter 'summation_small': %s.",
               summation_small);
        }
      }
    }
  }
}



NumericVector calculate_pdf(const int& Nrt, const int& Nres, const int& Na,
                            const int& Nv, const int& Nt0, const int& Nw,
                            const int& Nsv, const int& Nsig, const int& Nerr,
                            const int& Nmax,
                            const NumericVector& rt, const vector<int>& resp,
                            const NumericVector& a, const NumericVector& v,
                            const NumericVector& t0, const NumericVector& w,
                            const NumericVector& sv, const NumericVector& sigma,
                            const NumericVector& err,
                            const vector<bool>& invalid_input,
                            const int& max_terms_large,
                            const NumFunc& numf, const SumFunc& sumf,
                            const DenFunc& denf, const double& rt0)
{
  NumericVector out(Nmax);
  double t;
  if (Nsig == 1 && sigma[0] == 1) { // standard diffusion coefficient
    for (int i = 0; i < Nmax; i++) {
      t = rt[i % Nrt] - t0[i % Nt0]; // response time minus non-decision time
      // out-of-bounds handling, same priority as dnorm()
      if (invalid_input[i]) { // out-of-bounds model parameters
        out[i] = NAN;
        continue;
      }
      if (t <= 0) { // out-of-bounds response time
        out[i] = rt0;
        continue;
      }
      // sort response and calculate density
      if (resp[i % Nres] == 1) { // response is "lower" so use unchanged parameters
          out[i] = denf(t, a[i % Na], v[i % Nv], w[i % Nw], sv[i % Nsv],
                      err[i % Nerr], max_terms_large, numf, sumf);
      } else if (resp[i % Nres] == 2) { // response is "upper" so use alternate parameters
        out[i] = denf(t, a[i % Na], -v[i % Nv], 1 - w[i % Nw], sv[i % Nsv],
                      err[i % Nerr], max_terms_large, numf, sumf);
      } else { // out-of-bounds response, same output as o-o-b response time
        out[i] = rt0;
      }
    }
  } else { // non-standard diffusion coefficient
    for (int i = 0; i < Nmax; i++) {
      t = rt[i % Nrt] - t0[i % Nt0]; // response time minus non-decision time
      // out-of-bounds handling, same priority as dnorm()
      if (invalid_input[i]) { // out-of-bounds model parameters
        out[i] = NAN;
        continue;
      }
      if (t <= 0) { // out-of-bounds response time
        out[i] = rt0;
        continue;
      }
      // sort response and calculate density
      if (resp[i % Nres] == 1) { // response is "lower" so use unchanged parameters
        out[i] = denf(t, a[i % Na]/sigma[i % Nsig],
                      v[i % Nv]/sigma[i % Nsig], w[i % Nw],
                      sv[i % Nsv]/sigma[i % Nsig], err[i % Nerr],
                      max_terms_large, numf, sumf);
      } else if (resp[i % Nres] == 2) { // response is "upper" so use alternate parameters
        out[i] = denf(t, a[i % Na]/sigma[i % Nsig],
                      -v[i % Nv]/sigma[i % Nsig], 1 - w[i % Nw],
                      sv[i % Nsv]/sigma[i % Nsig], err[i % Nerr],
                      max_terms_large, numf, sumf);
      } else { // out-of-bounds response, same output as o-o-b response time
        out[i] = rt0;
      }
    }
  }

  return out;
}
