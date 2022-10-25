// Functions to check valid parameter values in DDM likelihood and gradient

#include "declarations.h"



vector<double> check_rt(const vector<double>& rt, int& Nrt) {
  Nrt = rt.size();
  vector<int> bad_idx;
  int bad_par = 0;

  for (int i = 0; i < Nrt; i++) {
    if (!(rt[i] > 0) || !isfinite(rt[i])) { // NaN, NA evaluate to !FALSE
      bad_par++;
      bad_idx.push_back(i);
    }
  }

  if (bad_par > 0) { // error handling
    string error_text = "fddm_fit error: response time is negative and/or infinite and/or NaN at the following indices: ";
    error_text = error_text.append(to_string(bad_idx[0]+1));
      for(int j = 1; j < bad_par; j++) {
        error_text = error_text.append(", ");
        error_text = error_text.append(to_string(bad_idx[j]+1));
      }
      error_text = error_text.append(".");
      stop(error_text);
  }

  return rt;
}


vector<double> convert_responses(const SEXP& response, const int& Nrt)
{
  vector<double> converted_responses;
  vector<int> bad_idx;
  int bad_par = 0;
  int type = TYPEOF(response);
  int Nres;
  
  if (type == 13 || type == 14) { // IntegerVector (including factors) or NumericVector
    converted_responses = Rcpp::as<vector<double> >(response);
    Nres = converted_responses.size();
    if (Nres < 1) {
      stop("fddm_fit error: function parameter 'response' is empty.");
    } else if (Nres == 1) { // convert response and duplicate until length Nrt
      converted_responses.resize(Nrt);
      double cr0 = converted_responses[0];
      if (!(cr0 == 1) && !(cr0 == 2)) { // NaN, NA evaluate to !FALSE; factors sort automatically
        bad_par++;
        bad_idx.push_back(0);
      }
      for (int i = 0; i < Nrt; i++) {
        converted_responses[i] = cr0;
      }
    } else if (Nres == Nrt) { // convert and copy each response
      converted_responses.resize(Nrt);
      double cri;
      for (int i = 0; i < Nrt; i++) {
        cri = converted_responses[i];
        if (!(cri == 1) && !(cri == 2)) { // NaN, NA evaluate to !FALSE
          bad_par++;
          bad_idx.push_back(i);
        }
      }
    } else {
      stop("fddm_fit error: number of responses is not 1 and does not match the number of response times");
    }
  } else if (type == 16) { // StringVector (contains at least one string)
    vector<string> temp = Rcpp::as<vector<string> >(response);
    Nres = temp.size();
    if (Nres < 1) {
      stop("fddm_fit error: function parameter 'response' is empty.");
    } else if (Nres == 1) {
      converted_responses.resize(Nrt);
      char temp0 = temp[0][0];
      if (temp0 == 'l' || temp0 == 'L') { // lower
        for (int i = 0; i < Nrt; i++) {
          converted_responses[i] = 1;
        }
      } else if (temp0 == 'u' || temp0 == 'U') { // upper
        for (int i = 0; i < Nrt; i++) {
          converted_responses[i] = 2;
        }
      } else {
        bad_par++;
        bad_idx.push_back(0);
      }
    } else if (Nres == Nrt) {
      converted_responses.resize(Nrt);
      for (int i = 0; i < Nrt; i++) {
        char tempi = temp[i][0];
        if (tempi == 'l' || tempi == 'L') {
          converted_responses[i] = 1;
        } else if (tempi == 'u' || tempi == 'U') {
          converted_responses[i] = 2;
        } else {
          bad_par++;
          bad_idx.push_back(i);
        }
      }
    } else {
      stop("fddm_fit error: number of responses is not 1 and does not match the number of response times");
    }
  } else if (type == 10) { // LogicalVector (boolean values)
    converted_responses = Rcpp::as<vector<double> >(response);
    Nres = converted_responses.size();
    if (Nres < 1) {
      stop("fddm_fit error: function parameter 'response' is empty.");
    } else if (Nres == 1) { // convert response and duplicate until length Nrt
      converted_responses.resize(Nrt);
      double cr0 = converted_responses[0] + 1;
      if (!(cr0 == 1) && !(cr0 == 2)) { // NaN, NA evaluate to !FALSE; factors sort automatically
        bad_par++;
        bad_idx.push_back(0);
      }
      for (int i = 0; i < Nrt; i++) {
        converted_responses[i] = cr0;
      }
    } else if (Nres == Nrt) { // convert and copy each response
      converted_responses.resize(Nrt);
      double cri;
      for (int i = 0; i < Nrt; i++) {
        cri = converted_responses[i] + 1;
        if (!(cri == 1) && !(cri == 2)) { // NaN, NA evaluate to !FALSE
          bad_par++;
          bad_idx.push_back(i);
        }
      }
    } else {
      stop("fddm_fit error: number of responses is not 1 and does not match the number of response times");
    }
  } else {
    stop("fddm_fit error: type of function parameter 'response' is not one of: integer, double, factor, string (character), or boolean (logical).");
  }
  
  if (bad_par > 0) { // error handling
    string error_text = "fddm_fit error: function parameter 'response' was input as a vector of ";
    
    if (type == 13) { // IntegerVector
      if (Rf_isFactor(response) == 1) { // factor
        error_text = error_text.append("factors, and a value that does not match either the first or second level");
      } else { // NOT factor
        error_text = error_text.append("integers, and a value other than 1 or 2");
      }
    } else if (type == 14) {
      error_text = error_text.append("doubles (truncated to integers), and a value other than 1 or 2");
    } else if (type == 16) {
      error_text = error_text.append("strings (characters), and a value other than 'L' or 'U' (case insensitive)");
    } else if (type == 10) {
      error_text = error_text.append("booleans (logicals), and a value other than TRUE or FALSE");
    } else {
      error_text = error_text.append("some unknown type");
    }
    
    error_text = error_text.append(" was detected at ");
    
    if (bad_par == 1) {
      error_text = error_text.append("index ");
      error_text = error_text.append(to_string(bad_idx[0]+1));
      error_text = error_text.append(".");
      stop(error_text);
    } else {
      error_text = error_text.append("the following indices: ");
      error_text = error_text.append(to_string(bad_idx[0]+1));
      for(int j = 1; j < bad_par; j++) {
        error_text = error_text.append(", ");
        error_text = error_text.append(to_string(bad_idx[j]+1));
      }
      error_text = error_text.append(".");
      stop(error_text);
    }
  }

  return converted_responses;
}


void unpack_and_check_mod_mats(const vector<MatrixXd>& model_matrices,
                               MatrixXd& mm_v, MatrixXd& mm_a, MatrixXd& mm_t0,
                               MatrixXd& mm_w, MatrixXd& mm_sv,
                               VectorXd& v, VectorXd& a, VectorXd& t0,
                               VectorXd& w, VectorXd& sv,
                               vector<int>& form_len, const int& Nrt) {
  // could maybe do a list of matrices and unpack them by name instead of position
  if (model_matrices.size() != 5) {
    stop("fddm_fit error: 5 model matrices were not supplied.");
  }

  mm_v = model_matrices[0];
  form_len[0] = mm_v.cols();
  if (form_len[0] == 1) {
    if (mm_v.rows() == 1) { // the DDM parameter is constant (not a formula)
      form_len[0] = 0;
      if (!isfinite(mm_v(0,0))) { // check constant `v`
        stop("fddm_fit error: DDM parameter `v` interpreted as a constant and is non-positive, infinite, or a NaN.");
      }
      for (int i = 0; i < Nrt; i++) {
        v[i] = mm_v(0,0);
      }
    } else if (mm_v.rows() != Nrt) {
      stop("fddm_fit error: the model matrix for DDM parameter `v` does not match the rt data.");
    }
  } // else it's a formula (perhaps with only one feature), regular matrix mult

  mm_a = model_matrices[1];
  form_len[1] = mm_a.cols();
  if (form_len[1] == 1) {
    if (mm_a.rows() == 1) { // the DDM parameter is constant (not a formula)
      form_len[1] = 0;
      if (!(mm_a(0,0) > 0) || !isfinite(mm_a(0,0))) { // check constant `a`
        stop("fddm_fit error: DDM parameter `a` interpreted as a constant and is non-positive, infinite, or a NaN.");
      }
      for (int i = 0; i < Nrt; i++) {
        a[i] = mm_a(0,0);
      }
    } else if (mm_a.rows() != Nrt) {
      stop("fddm_fit error: the model matrix for DDM parameter `a` does not match the rt data.");
    }
  } // else it's a formula (perhaps with only one feature), regular matrix mult
  
  mm_t0 = model_matrices[2];
  form_len[2] = mm_t0.cols();
  if (form_len[2] == 1) {
    if (mm_t0.rows() == 1) { // the DDM parameter is constant (not a formula)
      form_len[2] = 0;
      if (!(mm_t0(0,0) >= 0) || !isfinite(mm_t0(0,0))) { // check constant `a`
        stop("fddm_fit error: DDM parameter `t0` interpreted as a constant and is non-positive, infinite, or a NaN.");
      }
      for (int i = 0; i < Nrt; i++) {
        t0[i] = mm_t0(0,0);
      }
    } else if (mm_t0.rows() != Nrt) {
      stop("fddm_fit error: the model matrix for DDM parameter `t0` does not match the rt data.");
    }
  } // else it's a formula (perhaps with only one feature), regular matrix mult
  
  mm_w = model_matrices[3];
  form_len[3] = mm_w.cols();
  if (form_len[3] == 1) {
    if (mm_w.rows() == 1) { // the DDM parameter is constant (not a formula)
      form_len[3] = 0;
      if (!(mm_w(0,0) > 0) || !(mm_w(0,0) < 1)) { // check constant `w`
        stop("fddm_fit error: DDM parameter `w` interpreted as a constant and is non-positive, infinite, or a NaN.");
      }
      for (int i = 0; i < Nrt; i++) {
        w[i] = mm_w(0,0);
      }
    } else if (mm_w.rows() != Nrt) {
      stop("fddm_fit error: the model matrix for DDM parameter `w` does not match the rt data.");
    }
  } // else it's a formula (perhaps with only one feature), regular matrix mult
  
  mm_sv = model_matrices[4];
  form_len[4] = mm_sv.cols();
  if (form_len[4] == 1) {
    if (mm_sv.rows() == 1) { // the DDM parameter is constant (not a formula)
      form_len[4] = 0;
      if (!(mm_sv(0,0) >= 0) || !isfinite(mm_sv(0,0))) { // check constant `sv`
        stop("fddm_fit error: DDM parameter `sv` interpreted as a constant and is non-positive, infinite, or a NaN.");
      }
      for (int i = 0; i < Nrt; i++) {
        sv[i] = mm_sv(0,0);
      }
    } else if (mm_sv.rows() != Nrt) {
      stop("fddm_fit error: the model matrix for DDM parameter `sv` does not match the rt data.");
    }
  } // else it's a formula (perhaps with only one feature), regular matrix mult
}


double check_err_tol(const double& err_tol) {
  if (err_tol > 0 && isfinite(err_tol)) {
    if (err_tol > ERR_TOL_THRESH) {
      return err_tol;
    } else {
      warning("fddm_fit warning: function parameter 'err_tol' is positive, but it is too small; 'err_tol' has been set to %e.", ERR_TOL_THRESH);
      return ERR_TOL_THRESH;
    }
  } else {
    stop("fddm_fit error: function parameter 'err_tol' is non-positive and/or infinite: %e.", err_tol);
  }
}


double check_switch_thresh(const double& switch_thresh) {
  if (isnan(switch_thresh)) {
    stop("fddm_fit error: function parameter 'switch_thresh' is a NaN: %f.", switch_thresh);
  }
  return switch_thresh;
}


bool invalid_parameters(const VectorXd& v, const VectorXd& a,
                        const VectorXd& t0, const VectorXd& w,
                        const VectorXd& sv, const int& Nrt,
                        const vector<int>& form_len)
{
  // note: NaN, NA evaluate to FALSE and then get negated
  // if a parameter is constant, it was checked during construction

  if (form_len[0] > 0) {
    for (int i = 0; i < Nrt; i++) {
      if (!isfinite(v[i])) {
        return 1;
      }
    }
  }
  if (form_len[1] > 0) {
    for (int i = 0; i < Nrt; i++) {
      if (!(a[i] > 0) || !isfinite(a[i])) {
        return 1;
      }
    }
  }
  if (form_len[2] > 0) {
    for (int i = 0; i < Nrt; i++) {
      if (!(t0[i] >= 0) || !isfinite(t0[i])) {
        return 1;
      }
    }
  }
  if (form_len[3] > 0) {
    for (int i = 0; i < Nrt; i++) {
      if (!(w[i] > 0) || !(w[i] < 1)) {
        return 1;
      }
    }
  }
  if (form_len[4] > 0) {
    for (int i = 0; i < Nrt; i++) {
      if (!(sv[i] >= 0) || !isfinite(sv[i])) {
        return 1;
      }
    }
  }

  return 0;
}
