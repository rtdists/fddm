// Functions to check valid parameter values in dfddm and pfddm
// included from src/fitting_helper_functions/class_methods.h

using std::vector;
using std::string;
using std::to_string;
using std::isfinite;
using std::isnan;
using Rcpp::stop;
using Rcpp::warning;

// static const double ERR_TOL_THRESH = 1e-300; // near minimum value of a double


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
    string error_text = "ddm fit error: response time is negative and/or infinite and/or NaN at the following indices: ";
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
      stop("ddm fit error: function parameter 'response' is empty.");
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
      stop("ddm fit error: number of responses is not 1 and does not match the number of response times");
    }
  } else if (type == 16) { // StringVector (contains at least one string)
    vector<string> temp = Rcpp::as<vector<string> >(response);
    Nres = temp.size();
    if (Nres < 1) {
      stop("ddm fit error: function parameter 'response' is empty.");
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
      stop("ddm fit error: number of responses is not 1 and does not match the number of response times");
    }
  } else if (type == 10) { // LogicalVector (boolean values)
    converted_responses = Rcpp::as<vector<double> >(response);
    Nres = converted_responses.size();
    if (Nres < 1) {
      stop("ddm fit error: function parameter 'response' is empty.");
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
      stop("ddm fit error: number of responses is not 1 and does not match the number of response times");
    }
  } else {
    stop("ddm fit error: type of function parameter 'response' is not one of: integer, double, factor, string (character), or boolean (logical).");
  }
  
  if (bad_par > 0) { // error handling
    string error_text = "ddm fit error: function parameter 'response' was input as a vector of ";
    
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


double check_err_tol(const double& err_tol) {
  if (err_tol > 0 && isfinite(err_tol)) {
    if (err_tol > ERR_TOL_THRESH) {
      return err_tol;
    } else {
      warning("ddm fit warning: function parameter 'err_tol' is positive, but it is too small; 'err_tol' has been set to %e.", ERR_TOL_THRESH);
      return ERR_TOL_THRESH;
    }
  } else {
    stop("ddm fit error: function parameter 'err_tol' is non-positive and/or infinite: %e.", err_tol);
  }
}


double check_switch_thresh(const double& switch_thresh) {
  if (isnan(switch_thresh)) {
    stop("ddm fit error: function parameter 'switch_thresh' is a NaN: %f.", switch_thresh);
  }
  return switch_thresh;
}


bool invalid_parameters(vector<double>& parameters,
                        const vector<double>& temp_params)
{
  // unpack parameters, and store them in the class member 'parameters'
  parameters[0] = temp_params[0]; // a
  parameters[1] = temp_params[1]; // v
  parameters[2] = temp_params[2]; // t0
  parameters[3] = temp_params[3]; // w
  parameters[4] = temp_params[4]; // sv

  // note: NaN, NA evaluate to FALSE and then get negated
  if (!(parameters[0] > 0) || !isfinite(parameters[0])) {
    return 1;
  }
  
  if (!isfinite(parameters[1])) {
    return 1;
  }
  
  if (!(parameters[2] >= 0) || !isfinite(parameters[2])) {
    return 1;
  }
  
  if (!(parameters[3] > 0) || !(parameters[3] < 1)) {
    return 1;
  }

  if (!(parameters[4] > 0) || !isfinite(parameters[4])) {
    return 1;
  }
  
  return 0;
}
