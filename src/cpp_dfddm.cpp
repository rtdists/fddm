// Density function for the Ratcliff Diffusion Decision Model (DDM) PDF

#include <Rcpp.h>
// [[Rcpp::plugins(cpp11)]]
#define _USE_MATH_DEFINES
#include <math.h>
#include "funcs.h"




// General density function
// [[Rcpp::export]]
Rcpp::NumericVector cpp_dfddm(const Rcpp::NumericVector& rt,
                              const Rcpp::LogicalVector& response,
                              const double& a, const double& v,
                              const double& t0, const double& w,
                              const double& sv, const bool& log_prob,
                              const std::string& n_terms_small,
                              const std::string& summation_small,
                              const std::string& scale,
                              const double& eps)
{
  // switch-case through all of the options, maybe have a flag for the default?
  if (scale == "small") {
    if (summation_small == "2017") {
      if (n_terms_small == "Foster") {
        return fs_eps_2017(rt, response, a, v, t0, w, sv, log_prob, eps);
      } else if (n_terms_small == "Navarro") {
        return fs_Nav_2017(rt, response, a, v, t0, w, sv, log_prob, eps);
      } else if (n_terms_small == "Kesselmeier") {
        return fs_BGK_2017(rt, response, a, v, t0, w, sv, log_prob, eps);
      } else {
        Rcpp::Rcerr << "error: invalid n_terms_small" << std::endl;
        return NAN;
      }
    } else if (summation_small == "2014") {
      if (n_terms_small == "Foster") {
        return fs_eps_2014(rt, response, a, v, t0, w, sv, log_prob, eps);
      } else if (n_terms_small == "Navarro") {
        return fs_Nav_2014(rt, response, a, v, t0, w, sv, log_prob, eps);
      } else if (n_terms_small == "Kesselmeier") {
        return fs_BGK_2014(rt, response, a, v, t0, w, sv, log_prob, eps);
      } else {
        Rcpp::Rcerr << "error: invalid n_terms_small" << std::endl;
        return NAN;
      }
    } else {
      Rcpp::Rcerr << "error: invalid summation_small" << std::endl;
      return NAN;
    }
  } else if (scale == "large") {
    return fl_Nav_2009(rt, response, a, v, t0, w, log_prob, eps);
  } else if (scale == "both") {
    if (summation_small == "2017") {
      if (n_terms_small == "Navarro") {
        return fb_Nav_Nav_2017(rt, response, a, v, t0, w, sv, log_prob, eps);
      } else if (n_terms_small == "Kesselmeier") {
        return fb_BGK_Nav_2017(rt, response, a, v, t0, w, sv, log_prob, eps);
      } else {
        Rcpp::Rcerr << "error: invalid n_terms_small" << std::endl;
        return NAN;
      }
    } else if (summation_small == "2014") {
      if (n_terms_small == "Navarro") {
        return fb_Nav_Nav_2014(rt, response, a, v, t0, w, sv, log_prob, eps);
      } else if (n_terms_small == "Kesselmeier") {
        return fb_BGK_Nav_2014(rt, response, a, v, t0, w, sv, log_prob, eps);
      } else {
        Rcpp::Rcerr << "error: invalid n_terms_small" << std::endl;
        return NAN;
      }
    } else {
      Rcpp::Rcerr << "error: invalid summation_small" << std::endl;
      return NAN;
    }
  } else {
    Rcpp::Rcerr << "error: invalid scale" << std::endl;
    return NAN;
  }
}
