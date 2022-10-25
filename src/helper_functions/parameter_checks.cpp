#include "parameter_checks.h"


void convert_responses(const SEXP& response, int& Nres, int& Nmax,
                       vector<double>& out, const double& rt0, bool& valid)
{
  vector<int> bad_idx;
  int bad_par = 0;
  int type = TYPEOF(response);
  
  if (type == 13 || type == 14) { // IntegerVector (including factors) or NumericVector
    out = Rcpp::as<vector<double> >(response);
    Nres = out.size();
    if (Nres < 1) {
      warning("dfddm warning: function parameter 'response' is empty; empty vector returned.");
      valid = 0;
    } else {
      Nmax = max(Nmax, Nres);
      out.resize(Nmax);
      for (int i = 0; i < Nres; i++) {
        if (out[i] == 1 || out[i] == 2) { // factors sort automatically
          for (int j = i; j < Nmax; j += Nres) {
            out[j] = out[i];
          }
        } else { // response = {Inf, -Inf} implies PDF = 0 (or log(0) )
          double naan = isnan(out[i]) ? out[i] : rt0;
          for (int j = i; j < Nmax; j += Nres) {
            out[j] = naan;
          }
          bad_par++;
          bad_idx.push_back(i);
        }
      }
    }
  } else if (type == 16) { // StringVector (contains at least one string)
    vector<string> temp = Rcpp::as<vector<string> >(response);
    Nres = temp.size();
    if (Nres < 1) {
      warning("dfddm warning: function parameter 'response' is empty; empty vector returned.");
      valid = 0;
    } else {
      Nmax = max(Nmax, Nres);
      out.resize(Nmax);
      for (int i = 0; i < Nres; i++) {
        if (temp[i][0] == 'l' || temp[i][0] == 'L') { // lower
          for (int j = i; j < Nmax; j += Nres) {
            out[j] = 1;
          }
        } else if (temp[i][0] == 'u' || temp[i][0] == 'U') { // upper
          for (int j = i; j < Nmax; j += Nres) {
            out[j] = 2;
          }
        } else {
          double naan = rt0;
          if (temp[i][0] == 'N') {
            int temp_length = temp[i].size();
            if (temp_length == 2 && temp[i][1] == 'A') {
              naan = NA_REAL;
            } else if (temp_length == 3 && temp[i][1] == 'a' && temp[i][2] == 'N') {
              naan = NAN;
            }
          }
          for (int j = i; j < Nmax; j += Nres) {
            out[j] = naan;
          }
          bad_par++;
          bad_idx.push_back(i);
        }
      }
    }
  } else if (type == 10) { // LogicalVector (boolean values)
    out = Rcpp::as<vector<double> >(response);
    Nres = out.size();
    if (Nres < 1) {
      warning("dfddm warning: function parameter 'response' is empty; empty vector returned.");
      valid = 0;
    } else {
      Nmax = max(Nmax, Nres);
      out.reserve(Nmax);
      for (int i = 0; i < Nres; i++) {
        if (out[i] == 0 || out[i] == 1) {
          for (int j = i; j < Nmax; j += Nres) {
            out[j] = out[i] + 1;
          }
        } else { // response = {Inf, -Inf} implies PDF = 0 (or log(0) )
          double naan = isnan(out[i]) ? out[i] : rt0;
          for (int j = i; j < Nmax; j += Nres) {
            out[j] = naan;
          }
          bad_par++;
          bad_idx.push_back(i);
        }
      }
    }
  } else {
    stop("dfddm error: type of function parameter 'response' is not one of: integer, double, factor, string (character), or boolean (logical).");
  }
  
  if (bad_par > 0) { // error handling
    string warning_text = "dfddm warning: function parameter 'response' was input as a vector of ";
    
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
      warning_text = warning_text.append("booleans (logicals), and a value other than TRUE or FALSE");
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
      warning_text = warning_text.append("; zeros and/or NaNs produced.");
      warning(warning_text);
    }
  }
}



bool parameter_check(const int& Nrt, int& Nres, const int& Nv, const int& Na,
                     const int& Nt0, const int& Nw, const int& Nsv,
                     const int& Nsig, const int& Nerr, int& Nmax,
                     const NumericVector& rt, const SEXP& response,
                     const NumericVector& v, const NumericVector& a,
                     const NumericVector& t0, const NumericVector& w,
                     const NumericVector& sv, const NumericVector& sigma,
                     const NumericVector& err,
                     vector<double>& out, const double& rt0)
{
  bool valid = 1;
  
  // rt, invalid inputs will be handled in calculate_pdf()
  if (Nrt < 1) {
    warning("dfddm warning: function parameter 'rt' is empty; empty vector returned.");
    valid = 0;
  }
  
  // response, checks and converts responses
  convert_responses(response, Nres, Nmax, out, rt0, valid);
  
  // v
  if (Nv < 1) {
    warning("dfddm warning: model parameter 'v' is empty; empty vector returned.");
    valid = 0;
  } else {
    for (int i = 0; i < Nv; i++) {
      if (isfinite(v[i])) {
        continue;
      } else { // NaN, NA, Inf, -Inf are not finite
        double naan = isnan(v[i]) ? v[i] : rt0; // v = {Inf, -Inf} implies PDF = 0 (or log(0) )
        for (int j = i; j < Nmax; j += Nv) {
          out[j] = naan;
        }
      }
    }
  }
  
  // a
  if (Na < 1) {
    warning("dfddm warning: model parameter 'a' is empty; empty vector returned.");
    valid = 0;
  } else {
    for (int i = 0; i < Na; i++) {
      if (a[i] > 0) {
        if (isfinite(a[i])) {
          continue;
        } else { // a = Inf implies PDF = 0 (or log(0) )
          for (int j = i; j < Nmax; j += Na) {
            out[j] = rt0;
          }
        }
      } else { // {NaN, NA} evaluate to FALSE
        double naan = isnan(a[i]) ? a[i] : NAN;
        for (int j = i; j < Nmax; j += Na) {
          out[j] = naan;
        }
      }
    }
  }
  
  // t0
  if (Nt0 < 1) {
    warning("dfddm warning: model parameter 't0' is empty; empty vector returned.");
    valid = 0;
  } else {
    for (int i = 0; i < Nt0; i++) {
      if (t0[i] >= 0) {
        if (isfinite(t0[i])) { // this could also be handled in calculate_pdf()
          continue;
        } else { // t0 = Inf implies rt - t0 < 0 implies PDF = 0 (or log(0) )
          for (int j = i; j < Nmax; j += Nt0) {
            out[j] = rt0;
          }
        }
      } else { // {NaN, NA} evaluate to FALSE
        double naan = isnan(t0[i]) ? t0[i] : NAN;
        for (int j = i; j < Nmax; j += Nt0) {
          out[j] = naan;
        }
      }
    }
  }
  
  // w
  if (Nw < 1) {
    warning("dfddm warning: model parameter 'w' is empty; empty vector returned.");
    valid = 0;
  } else {
    for (int i = 0; i < Nw; i++) {
      if (w[i] > 0 && w[i] < 1) {
        continue;
      } else { // {NaN, NA} evaluate to FALSE
        double naan = isnan(w[i]) ? w[i] : NAN;
        for (int j = i; j < Nmax; j += Nw) {
          out[j] = naan;
        }
      }
    }
  }
  
  // sv
  if (Nsv < 1) {
    warning("dfddm warning: model parameter 'sv' is empty; empty vector returned.");
    valid = 0;
  } else {
    for (int i = 0; i < Nsv; i++) {
      if (sv[i] >= 0) {
        if (isfinite(sv[i])) {
          continue;
        } else { // sv = Inf implies PDF = 0 (or log(0) )
          for (int j = i; j < Nmax; j += Nsv) {
            out[j] = rt0;
          }
        }
      } else { // {NaN, NA} evaluate to FALSE
        double naan = isnan(sv[i]) ? sv[i] : NAN;
        for (int j = i; j < Nmax; j += Nsv) {
          out[j] = naan;
        }
      }
    }
  }
  
  // sigma
  if (Nsig < 1) {
    warning("dfddm warning: model parameter 'sigma' is empty; empty vector returned.");
    valid = 0;
  } else {
    for (int i = 0; i < Nsig; i++) {
      if (sigma[i] > 0 && isfinite(sigma[i])) { // sigma = Inf causes problems
        continue;
      } else { // {NaN, NA} evaluate to FALSE
        double naan = isnan(sigma[i]) ? sigma[i] : NAN;
        for (int j = i; j < Nmax; j += Nsig) {
          out[j] = naan;
        }
      }
    }
  }
  
  // err_tol
  if (Nerr < 1) {
    stop("dfddm error: function parameter 'err_tol' is empty.");
    valid = 0;
  } else {
    int bad_par = 0;
    vector<int> bad_idx;
    for (int i = 0; i < Nerr; i++) {
      if (err[i] > 0) {
        if (isfinite(err[i])) {
          continue;
        } else {
          bad_par++;
          bad_idx.push_back(i);
        }
      } else { // {NaN, NA} evaluate to FALSE
        bad_par++;
        bad_idx.push_back(i);
      }
    }
    if (bad_par > 0) {
      if (bad_par == 1) {
        stop("dfddm error: function parameter 'err_tol' <= 0 or is not finite at index %i.",
             bad_idx[0]+1);
      } else {
        std::string error_text = "dfddm error: function parameter 'err_tol' <= 0 or is not finite at the following indices:\n";
        error_text = error_text.append(to_string(bad_idx[0]+1));
        for(int j = 1; j < bad_par; j++) {
          error_text = error_text.append(", ");
          error_text = error_text.append(to_string(bad_idx[j]+1));
        }
        error_text = error_text.append(".");
        stop(error_text);
      }
    }
  }
  
  return valid;
}
