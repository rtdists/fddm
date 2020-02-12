// Density function for the Ratcliff Diffusion Decision Model (DDM) PDF

#include "funcs.h"

using std::endl;
using Rcpp::Rcerr;




// General density function
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
  // if (min(rt.length(),
  //         response.length(),
  //         a.length(),
  //         v.length(),
  //         t0.length(),
  //         w.length(),
  //         sv.length(),
  //         log_prob.length(),
  //         eps.length()) < 1) {
  //   Rcerr << "error: one of the inputs does not exist" << endl;
  //   return NAN;
  // }
  
  int Nmax = max({rt.length(),
                 response.length(),
                 a.length(),
                 v.length(),
                 t0.length(),
                 w.length(),
                 sv.length(),
                 log_prob.length(),
                 eps.length()});
  return Nmax;
  
  // if (!log_prob) { // non-log version
  //   if (scale == "small") {
  //     if (summation_small == "2017") {
  //       if (n_terms_small == "Foster") {
  //         return fs_eps_2017(rt, response, a, v, t0, w, sv, eps);
  //       } else if (n_terms_small == "Navarro") {
  //         return fs_Nav_2017(rt, response, a, v, t0, w, sv, eps);
  //       } else if (n_terms_small == "Kesselmeier") {
  //         return fs_BGK_2017(rt, response, a, v, t0, w, sv, eps);
  //       } else {
  //         Rcerr << "error: invalid n_terms_small" << endl;
  //         return NAN;
  //       }
  //     } else if (summation_small == "2014") {
  //       if (n_terms_small == "Foster") {
  //         return fs_eps_2014(rt, response, a, v, t0, w, sv, eps);
  //       } else if (n_terms_small == "Navarro") {
  //         return fs_Nav_2014(rt, response, a, v, t0, w, sv, eps);
  //       } else if (n_terms_small == "Kesselmeier") {
  //         return fs_BGK_2014(rt, response, a, v, t0, w, sv, eps);
  //       } else {
  //         Rcerr << "error: invalid n_terms_small" << endl;
  //         return NAN;
  //       }
  //     } else {
  //       Rcerr << "error: invalid summation_small" << endl;
  //       return NAN;
  //     }
  //   } else if (scale == "large") {
  //     return fl_Nav_2009(rt, response, a, v, t0, w, eps);
  //   } else if (scale == "both") {
  //     if (summation_small == "2017") {
  //       if (n_terms_small == "Navarro") {
  //         return fb_Nav_Nav_2017(rt, response, a, v, t0, w, sv, eps);
  //       } else if (n_terms_small == "Kesselmeier") {
  //         return fb_BGK_Nav_2017(rt, response, a, v, t0, w, sv, eps);
  //       } else {
  //         Rcerr << "error: invalid n_terms_small" << endl;
  //         return NAN;
  //       }
  //     } else if (summation_small == "2014") {
  //       if (n_terms_small == "Navarro") {
  //         return fb_Nav_Nav_2014(rt, response, a, v, t0, w, sv, eps);
  //       } else if (n_terms_small == "Kesselmeier") {
  //         return fb_BGK_Nav_2014(rt, response, a, v, t0, w, sv, eps);
  //       } else {
  //         Rcerr << "error: invalid n_terms_small" << endl;
  //         return NAN;
  //       }
  //     } else {
  //       Rcerr << "error: invalid summation_small" << endl;
  //       return NAN;
  //     }
  //   } else {
  //     Rcerr << "error: invalid scale" << endl;
  //     return NAN;
  //   }
  // } else { // log version
  //   if (scale == "small") {
  //     if (summation_small == "2017") {
  //       if (n_terms_small == "Foster") {
  //         return fs_eps_2017_log(rt, response, a, v, t0, w, sv, eps);
  //       } else if (n_terms_small == "Navarro") {
  //         return fs_Nav_2017_log(rt, response, a, v, t0, w, sv, eps);
  //       } else if (n_terms_small == "Kesselmeier") {
  //         return fs_BGK_2017_log(rt, response, a, v, t0, w, sv, eps);
  //       } else {
  //         Rcerr << "error: invalid n_terms_small" << endl;
  //         return NAN;
  //       }
  //     } else if (summation_small == "2014") {
  //       if (n_terms_small == "Foster") {
  //         return fs_eps_2014_log(rt, response, a, v, t0, w, sv, eps);
  //       } else if (n_terms_small == "Navarro") {
  //         return fs_Nav_2014_log(rt, response, a, v, t0, w, sv, eps);
  //       } else if (n_terms_small == "Kesselmeier") {
  //         return fs_BGK_2014_log(rt, response, a, v, t0, w, sv, eps);
  //       } else {
  //         Rcerr << "error: invalid n_terms_small" << endl;
  //         return NAN;
  //       }
  //     } else {
  //       Rcerr << "error: invalid summation_small" << endl;
  //       return NAN;
  //     }
  //   } else if (scale == "large") {
  //     return fl_Nav_2009_log(rt, response, a, v, t0, w, eps);
  //   } else if (scale == "both") {
  //     if (summation_small == "2017") {
  //       if (n_terms_small == "Navarro") {
  //         return fb_Nav_Nav_2017_log(rt, response, a, v, t0, w, sv, eps);
  //       } else if (n_terms_small == "Kesselmeier") {
  //         return fb_BGK_Nav_2017_log(rt, response, a, v, t0, w, sv, eps);
  //       } else {
  //         Rcerr << "error: invalid n_terms_small" << endl;
  //         return NAN;
  //       }
  //     } else if (summation_small == "2014") {
  //       if (n_terms_small == "Navarro") {
  //         return fb_Nav_Nav_2014_log(rt, response, a, v, t0, w, sv, eps);
  //       } else if (n_terms_small == "Kesselmeier") {
  //         return fb_BGK_Nav_2014_log(rt, response, a, v, t0, w, sv, eps);
  //       } else {
  //         Rcerr << "error: invalid n_terms_small" << endl;
  //         return NAN;
  //       }
  //     } else {
  //       Rcerr << "error: invalid summation_small" << endl;
  //       return NAN;
  //     }
  //   } else {
  //     Rcerr << "error: invalid scale" << endl;
  //     return NAN;
  //   }
  // }
}
