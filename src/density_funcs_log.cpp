// Functions to evaluate the log of the DDM PDF for specific criteria

#include <Rcpp.h>
// [[Rcpp::plugins(cpp11)]]
#define _USE_MATH_DEFINES
#include <math.h>
#include "funcs.h"




//////////                                                           //////////
///////////////////////////////// Small Time //////////////////////////////////
//////////                                                           //////////

// Use term < eps BGK2017 style sum approximation, with minimum terms, log
Rcpp::NumericVector fs_eps_2017_log(const Rcpp::NumericVector& rt,
                                    Rcpp::LogicalVector response,
                                    const double& a, const double& v,
                                    const double& t0, const double& w,
                                    const double& sv, const double& eps)
{
  int n = rt.length(); // get number of response times

  if (response.length() != n) { // create valid binary response vector
    bool first_resp = response[0];
    Rcpp::LogicalVector resp(n, first_resp);
    response = resp;
  }

  Rcpp::NumericVector out(n);
  double t, mult;
  double vprime = -v;
  double wprime = 1 - w;
  double log_a = log(a);

  if (sv < SV_THRESH) { // set sv=0 and use constant drift rate method
    for (int i = 0; i < n; i++) { // iterate through all values in rt and out
      t = rt[i] - t0; // subtract non-decisison time from the response time
      if (t <= 0) {
        out[i] = 0; // if t=0, return 0 instead of NaN
        continue;
      }

      double log_t3 = 1.5 * log(t);
      if (response[i] == 1) { // if response is "upper" use alternate parameters
        mult = log_a - LOG_2PI_2 - log_t3 - vprime * a * wprime
             - vprime * vprime * t / 2;
        out[i] = mult + log(small_sum_eps_17(t, a, wprime, eps/exp(mult)));
      } else { // else response is "lower"
        mult = log_a - LOG_2PI_2 - log_t3 -v * a * w - v * v * t / 2;
        out[i] = mult + log(small_sum_eps_17(t, a, w, eps/exp(mult)));
      }
    }
  } else { // use variable drift rate (changes mult)
    for (int i = 0; i < n; i++) { // iterate through all values in rt and out
      t = rt[i] - t0; // subtract non-decisison time from the response time
      if (t <= 0) {
        out[i] = 0; // if t=0, return 0 instead of NaN
        continue;
      }

      double log_t3 = 1.5 * log(t);
      double log_svt = 0.5 * log(1 + sv * t);
      if (response[i] == 1) { // if response is "upper" use alternate parameters
        mult = log_a - log_t3 - log_svt + (sv * a * a * wprime * wprime
             - 2 * vprime * a * wprime - vprime * vprime * t) / (2 + 2 * sv * t);
        out[i] = mult + log(small_sum_eps_17(t, a, wprime, eps/exp(mult)));
      } else { // else response is "lower"
        mult = log_a - log_t3 - log_svt + (sv * a * a * w * w
             - 2 * v * a * w - v * v * t) / (2 + 2 * sv * t);
        out[i] = mult + log(small_sum_eps_17(t, a, w, eps/exp(mult)));
      }
    }
  }
  return out;
}


// Use term < eps BGK2014 style sum approximation, with minimum terms, log
Rcpp::NumericVector fs_eps_2014_log(const Rcpp::NumericVector& rt,
                                    Rcpp::LogicalVector response,
                                    const double& a, const double& v,
                                    const double& t0, const double& w,
                                    const double& sv, const double& eps)
{
  int n = rt.length(); // get number of response times

  if (response.length() != n) { // create valid binary response vector
    bool first_resp = response[0];
    Rcpp::LogicalVector resp(n, first_resp);
    response = resp;
  }

  Rcpp::NumericVector out(n);
  double t, mult;
  double vprime = -v;
  double wprime = 1 - w;
  double log_a = log(a);

  if (sv < SV_THRESH) { // set sv=0 and use constant drift rate method
    for (int i = 0; i < n; i++) { // iterate through all values in rt and out
      t = rt[i] - t0; // subtract non-decisison time from the response time
      if (t <= 0) {
        out[i] = 0; // if t=0, return 0 instead of NaN
        continue;
      }

      double log_t3 = 1.5 * log(t);
      if (response[i] == 1) { // if response is "upper" use alternate parameters
        mult = log_a - LOG_2PI_2 - log_t3 - vprime * a * wprime
             - vprime * vprime * t / 2;
        out[i] = mult + log(small_sum_eps_14(t, a, wprime, eps/exp(mult)));
      } else { // else response is "lower"
        mult = log_a - LOG_2PI_2 - log_t3 -v * a * w - v * v * t / 2;
        out[i] = mult + log(small_sum_eps_14(t, a, w, eps/exp(mult)));
      }
    }
  } else { // use variable drift rate (changes mult)
    for (int i = 0; i < n; i++) { // iterate through all values in rt and out
      t = rt[i] - t0; // subtract non-decisison time from the response time
      if (t <= 0) {
        out[i] = 0; // if t=0, return 0 instead of NaN
        continue;
      }

      double log_t3 = 1.5 * log(t);
      double log_svt = 0.5 * log(1 + sv * t);
      if (response[i] == 1) { // if response is "upper" use alternate parameters
        mult = log_a - log_t3 - log_svt + (sv * a * a * wprime * wprime
             - 2 * vprime * a * wprime - vprime * vprime * t) / (2 + 2 * sv * t);
        out[i] = mult + log(small_sum_eps_14(t, a, wprime, eps/exp(mult)));
      } else { // else response is "lower"
        mult = log_a - log_t3 - log_svt + (sv * a * a * w * w
             - 2 * v * a * w - v * v * t) / (2 + 2 * sv * t);
        out[i] = mult + log(small_sum_eps_14(t, a, w, eps/exp(mult)));
      }
    }
  }
  return out;
}


// Use Navarro2009 number of terms for 2017 style sum approximation, log
Rcpp::NumericVector fs_Nav_2017_log(const Rcpp::NumericVector& rt,
                                    Rcpp::LogicalVector response,
                                    const double& a, const double& v,
                                    const double& t0, const double& w,
                                    const double& sv, const double& eps)
{
  int n = rt.length(); // get number of response times

  if (response.length() != n) { // create valid binary response vector
    bool first_resp = response[0];
    Rcpp::LogicalVector resp(n, first_resp);
    response = resp;
  }

  Rcpp::NumericVector out(n);
  double t, mult;
  int ks;
  double vprime = -v;
  double wprime = 1 - w;
  double log_a = log(a);

  if (sv < SV_THRESH) { // set sv=0 and use constant drift rate method
    for (int i = 0; i < n; i++) { // iterate through all values in rt and out
      t = rt[i] - t0; // subtract non-decisison time from the response time
      if (t <= 0) {
        out[i] = 0; // if t=0, return 0 instead of NaN
        continue;
      }

      ks = ks_Nav(t / (a * a), eps); // number of terms in sum approximation
      double log_t3 = 1.5 * log(t);
      if (response[i] == 1) { // if response is "upper" use alternate parameters
        mult = log_a - LOG_2PI_2 - log_t3 - vprime * a * wprime
             - vprime * vprime * t / 2;
        out[i] = mult + log(small_sum_2017(t, a, wprime, ks));
      } else { // else response is "lower"
        mult = log_a - LOG_2PI_2 - log_t3 -v * a * w - v * v * t / 2;
        out[i] = mult + log(small_sum_2017(t, a, w, ks));
      }
    }
  } else { // use variable drift rate (changes mult)
    for (int i = 0; i < n; i++) { // iterate through all values in rt and out
      t = rt[i] - t0; // subtract non-decisison time from the response time
      if (t <= 0) {
        out[i] = 0; // if t=0, return 0 instead of NaN
        continue;
      }

      ks = ks_Nav(t / (a * a), eps); // number of terms in sum approximation
      double log_t3 = 1.5 * log(t);
      double log_svt = 0.5 * log(1 + sv * t);
      if (response[i] == 1) { // if response is "upper" use alternate parameters
        mult = log_a - log_t3 - log_svt + (sv * a * a * wprime * wprime
             - 2 * vprime * a * wprime - vprime * vprime * t) / (2 + 2 * sv * t);
        out[i] = mult + log(small_sum_2017(t, a, wprime, ks));
      } else { // else response is "lower"
        mult = log_a - log_t3 - log_svt + (sv * a * a * w * w
             - 2 * v * a * w - v * v * t) / (2 + 2 * sv * t);
        out[i] = mult + log(small_sum_2017(t, a, w, ks));
      }
    }
  }
  return out;
}


// Use Navarro2009 number of terms for 2014 style sum approximation, log
Rcpp::NumericVector fs_Nav_2014_log(const Rcpp::NumericVector& rt,
                                    Rcpp::LogicalVector response,
                                    const double& a, const double& v,
                                    const double& t0, const double& w,
                                    const double& sv, const double& eps)
{
  int n = rt.length(); // get number of response times

  if (response.length() != n) { // create valid binary response vector
    bool first_resp = response[0];
    Rcpp::LogicalVector resp(n, first_resp);
    response = resp;
  }

  Rcpp::NumericVector out(n);
  double t, mult;
  int ks;
  double vprime = -v;
  double wprime = 1 - w;
  double log_a = log(a);

  if (sv < SV_THRESH) { // set sv=0 and use constant drift rate method
    for (int i = 0; i < n; i++) { // iterate through all values in rt and out
      t = rt[i] - t0; // subtract non-decisison time from the response time
      if (t <= 0) {
        out[i] = 0; // if t=0, return 0 instead of NaN
        continue;
      }

      ks = ks_Nav(t / (a * a), eps); // number of terms in sum approximation
      double log_t3 = 1.5 * log(t);
      if (response[i] == 1) { // if response is "upper" use alternate parameters
        mult = log_a - LOG_2PI_2 - log_t3 - vprime * a * wprime
             - vprime * vprime * t / 2;
        out[i] = mult + log(small_sum_2014(t, a, wprime, ks));
      } else { // else response is "lower"
        mult = log_a - LOG_2PI_2 - log_t3 -v * a * w - v * v * t / 2;
        out[i] = mult + log(small_sum_2014(t, a, w, ks));
      }
    }
  } else { // use variable drift rate (changes mult)
    for (int i = 0; i < n; i++) { // iterate through all values in rt and out
      t = rt[i] - t0; // subtract non-decisison time from the response time
      if (t <= 0) {
        out[i] = 0; // if t=0, return 0 instead of NaN
        continue;
      }

      ks = ks_Nav(t / (a * a), eps); // number of terms in sum approximation
      double log_t3 = 1.5 * log(t);
      double log_svt = 0.5 * log(1 + sv * t);
      if (response[i] == 1) { // if response is "upper" use alternate parameters
        mult = log_a - log_t3 - log_svt + (sv * a * a * wprime * wprime
             - 2 * vprime * a * wprime - vprime * vprime * t) / (2 + 2 * sv * t);
        out[i] = mult + log(small_sum_2014(t, a, wprime, ks));
      } else { // else response is "lower"
        mult = log_a - log_t3 - log_svt + (sv * a * a * w * w
             - 2 * v * a * w - v * v * t) / (2 + 2 * sv * t);
        out[i] = mult + log(small_sum_2014(t, a, w, ks));
      }
    }
  }
  return out;
}


// Use BGK2014 number of terms for 2017 style sum approximation
Rcpp::NumericVector fs_BGK_2017_log(const Rcpp::NumericVector& rt,
                                    Rcpp::LogicalVector response,
                                    const double& a, const double& v,
                                    const double& t0, const double& w,
                                    const double& sv, const double& eps)
{
  int n = rt.length(); // get number of response times

  if (response.length() != n) { // create valid binary response vector
    bool first_resp = response[0];
    Rcpp::LogicalVector resp(n, first_resp);
    response = resp;
  }

  Rcpp::NumericVector out(n);
  double t, mult;
  int ks;
  double vprime = -v;
  double wprime = 1 - w;
  double log_a = log(a);

  if (sv < SV_THRESH) { // set sv=0 and use constant drift rate method
    for (int i = 0; i < n; i++) { // iterate through all values in rt and out
      t = rt[i] - t0; // subtract non-decisison time from the response time
      if (t <= 0) {
        out[i] = 0; // if t=0, return 0 instead of NaN
        continue;
      }

      double log_t3 = 1.5 * log(t);
      if (response[i] == 1) { // if response is "upper" use alternate parameters
        ks = ks_BGK(t / (a * a), wprime, eps); // number of terms in sum approximation
        mult = log_a - LOG_2PI_2 - log_t3 - vprime * a * wprime
             - vprime * vprime * t / 2;
        out[i] = mult + log(small_sum_2017(t, a, wprime, ks));
      } else { // else response is "lower"
        ks = ks_BGK(t / (a * a), w, eps); // number of terms in sum approximation
        mult = log_a - LOG_2PI_2 - log_t3 -v * a * w - v * v * t / 2;
        out[i] = mult + log(small_sum_2017(t, a, w, ks));
      }
    }
  } else { // use variable drift rate (changes mult)
    for (int i = 0; i < n; i++) { // iterate through all values in rt and out
      t = rt[i] - t0; // subtract non-decisison time from the response time
      if (t <= 0) {
        out[i] = 0; // if t=0, return 0 instead of NaN
        continue;
      }

      double log_t3 = 1.5 * log(t);
      double log_svt = 0.5 * log(1 + sv * t);
      if (response[i] == 1) { // if response is "upper" use alternate parameters
        ks = ks_BGK(t / (a * a), wprime, eps); // number of terms in sum approximation
        mult = log_a - log_t3 - log_svt + (sv * a * a * wprime * wprime
             - 2 * vprime * a * wprime - vprime * vprime * t) / (2 + 2 * sv * t);
        out[i] = mult + log(small_sum_2017(t, a, wprime, ks));
      } else { // else response is "lower"
        ks = ks_BGK(t / (a * a), w, eps); // number of terms in sum approximation
        mult = log_a - log_t3 - log_svt + (sv * a * a * w * w
             - 2 * v * a * w - v * v * t) / (2 + 2 * sv * t);
        out[i] = mult + log(small_sum_2017(t, a, w, ks));
      }
    }
  }
  return out;
}


// Use BGK2014 number of terms for 2017 style sum approximation
Rcpp::NumericVector fs_BGK_2014_log(const Rcpp::NumericVector& rt,
                                    Rcpp::LogicalVector response,
                                    const double& a, const double& v,
                                    const double& t0, const double& w,
                                    const double& sv, const double& eps)
{
  int n = rt.length(); // get number of response times

  if (response.length() != n) { // create valid binary response vector
    bool first_resp = response[0];
    Rcpp::LogicalVector resp(n, first_resp);
    response = resp;
  }

  Rcpp::NumericVector out(n);
  double t, mult;
  int ks;
  double vprime = -v;
  double wprime = 1 - w;
  double log_a = log(a);

  if (sv < SV_THRESH) { // set sv=0 and use constant drift rate method
    for (int i = 0; i < n; i++) { // iterate through all values in rt and out
      t = rt[i] - t0; // subtract non-decisison time from the response time
      if (t <= 0) {
        out[i] = 0; // if t=0, return 0 instead of NaN
        continue;
      }

      double log_t3 = 1.5 * log(t);
      if (response[i] == 1) { // if response is "upper" use alternate parameters
        ks = ks_BGK(t / (a * a), wprime, eps); // number of terms in sum approximation
        mult = log_a - LOG_2PI_2 - log_t3 - vprime * a * wprime
             - vprime * vprime * t / 2;
        out[i] = mult + log(small_sum_2014(t, a, wprime, ks));
      } else { // else response is "lower"
        ks = ks_BGK(t / (a * a), w, eps); // number of terms in sum approximation
        mult = log_a - LOG_2PI_2 - log_t3 -v * a * w - v * v * t / 2;
        out[i] = mult + log(small_sum_2014(t, a, w, ks));
      }
    }
  } else { // use variable drift rate (changes mult)
    for (int i = 0; i < n; i++) { // iterate through all values in rt and out
      t = rt[i] - t0; // subtract non-decisison time from the response time
      if (t <= 0) {
        out[i] = 0; // if t=0, return 0 instead of NaN
        continue;
      }

      double log_t3 = 1.5 * log(t);
      double log_svt = 0.5 * log(1 + sv * t);
      if (response[i] == 1) { // if response is "upper" use alternate parameters
        ks = ks_BGK(t / (a * a), wprime, eps); // number of terms in sum approximation
        mult = log_a - log_t3 - log_svt + (sv * a * a * wprime * wprime
             - 2 * vprime * a * wprime - vprime * vprime * t) / (2 + 2 * sv * t);
        out[i] = mult + log(small_sum_2014(t, a, wprime, ks));
      } else { // else response is "lower"
        ks = ks_BGK(t / (a * a), w, eps); // number of terms in sum approximation
        mult = log_a - log_t3 - log_svt + (sv * a * a * w * w
             - 2 * v * a * w - v * v * t) / (2 + 2 * sv * t);
        out[i] = mult + log(small_sum_2014(t, a, w, ks));
      }
    }
  }
  return out;
}



//////////                                                           //////////
///////////////////////////////// Large Time //////////////////////////////////
//////////                                                           //////////

// Use Navarro2009 number of terms for sum approximation
Rcpp::NumericVector fl_Nav_2009_log(const Rcpp::NumericVector& rt,
                                    Rcpp::LogicalVector response,
                                    const double& a, const double& v,
                                    const double& t0, const double& w,
                                    const double& eps)
{
  int n = rt.length(); // get number of response times

  if (response.length() != n) { // create valid binary response vector
    bool first_resp = response[0];
    Rcpp::LogicalVector resp(n, first_resp);
    response = resp;
  }

  Rcpp::NumericVector out(n);
  double t, mult;
  int kl;
  double vprime = -v;
  double wprime = 1 - w;
  double log_a2 = 2 * log(a);

  for (int i = 0; i < n; i++) { // iterate through all values in rt and out
    t = rt[i] - t0; // subtract non-decisison time from the response time
    if (t <= 0) {
      out[i] = 0; // if t=0, return 0 instead of NaN
      continue;
    }

    kl = kl_Nav(t / (a * a), eps); // number of terms in sum approximation
    if (response[i] == 1) { // if response is "upper" use alternate parameters
      mult = LOG_PI - log_a2 - vprime * a * wprime - vprime * vprime * t / 2;
      out[i] = mult + log(large_sum_Nav(t, a, wprime, kl));
    } else { // else response is "lower"
      mult =  LOG_PI - log_a2 - v * a * w - v * v * t / 2;
      out[i] = mult + log(large_sum_Nav(t, a, w, kl));
    }
  }
  return out;
}



//////////                                                           //////////
///////////////////////////////// Both Times //////////////////////////////////
//////////                                                           //////////

// ks = Navarro2009, 2017 style sum approximation, kl = Navarro2009
Rcpp::NumericVector fb_Nav_Nav_2017_log(const Rcpp::NumericVector& rt,
                                        Rcpp::LogicalVector response,
                                        const double& a, const double& v,
                                        const double& t0, const double& w,
                                        const double& sv, const double& eps)
{
  int n = rt.length(); // get number of response times

  if (response.length() != n) { // create valid binary response vector
    bool first_resp = response[0];
    Rcpp::LogicalVector resp(n, first_resp);
    response = resp;
  }

  Rcpp::NumericVector out(n);
  double t, mult;
  double vprime = -v;
  double wprime = 1 - w;

  if (sv < SV_THRESH) { // set sv=0 and use constant drift rate methods
    double log_a = log(a);
    double log_a2 = 2 * log_a;
    int ks, kl;
    for (int i = 0; i < n; i++) { // iterate through all values in rt and out
      t = rt[i] - t0; // subtract non-decisison time from the response time
      if (t <= 0) {
        out[i] = 0; // if t=0, return 0 instead of NaN
        continue;
      }

      ks = ks_Nav(t / (a * a), eps); // number of terms in small time sum approximation
      kl = kl_Nav(t / (a * a), eps); // number of terms in sum approximation
      double log_t3 = 1.5 * log(t);
      if (response[i] == 1) { // if response is "upper" use alternate parameters
        if (ks < kl) { // small time needs fewer terms than lorge time
          mult = log_a - LOG_2PI_2 - log_t3 - vprime * a * wprime
               - vprime * vprime * t / 2;
          out[i] = mult + log(small_sum_2017(t, a, wprime, ks));
        } else { // large time needs fewer terms
          mult = LOG_PI - log_a2 - vprime * a * wprime - vprime * vprime * t / 2;
          out[i] = mult + log(large_sum_Nav(t, a, wprime, kl));
        }
      } else { // else response is "lower"
        if (ks < kl) { // small time needs fewer terms than lorge time
          mult = log_a - LOG_2PI_2 - log_t3 -v * a * w - v * v * t / 2;
          out[i] = mult + log(small_sum_2017(t, a, w, ks));
        } else { // large time needs fewer terms
          mult =  LOG_PI - log_a2 - v * a * w - v * v * t / 2;
          out[i] = mult + log(large_sum_Nav(t, a, w, kl));
        }
      }
    }
  }
  else { // use variable drift rate (changes mult) and only use small time
    double log_a = log(a);
    int ks;
    for (int i = 0; i < n; i++) { // iterate through all values in rt and out
      t = rt[i] - t0; // subtract non-decisison time from the response time
      if (t <= 0) {
        out[i] = 0; // if t=0, return 0 instead of NaN
        continue;
      }

      ks = ks_Nav(t / (a * a), eps); // number of terms in sum approximation
      double log_t3 = 1.5 * log(t);
      double log_svt = 0.5 * log(1 + sv * t);
      if (response[i] == 1) { // if response is "upper" use alternate parameters
        mult = log_a - log_t3 - log_svt + (sv * a * a * wprime * wprime
             - 2 * vprime * a * wprime - vprime * vprime * t) / (2 + 2 * sv * t);
        out[i] = mult + log(small_sum_2017(t, a, wprime, ks));
      } else { // else response is "lower"
        mult = log_a - log_t3 - log_svt + (sv * a * a * w * w
             - 2 * v * a * w - v * v * t) / (2 + 2 * sv * t);
        out[i] = mult + log(small_sum_2017(t, a, w, ks));
      }
    }
  }
  return out;
}


// ks = Navarro2009, 2014 style sum approximation, kl = Navarro2009
Rcpp::NumericVector fb_Nav_Nav_2014_log(const Rcpp::NumericVector& rt,
                                        Rcpp::LogicalVector response,
                                        const double& a, const double& v,
                                        const double& t0, const double& w,
                                        const double& sv, const double& eps)
{
  int n = rt.length(); // get number of response times

  if (response.length() != n) { // create valid binary response vector
    bool first_resp = response[0];
    Rcpp::LogicalVector resp(n, first_resp);
    response = resp;
  }

  Rcpp::NumericVector out(n);
  double t, mult;
  double vprime = -v;
  double wprime = 1 - w;

  if (sv < SV_THRESH) { // set sv=0 and use constant drift rate methods
    double log_a = log(a);
    double log_a2 = 2 * log_a;
    int ks, kl;
    for (int i = 0; i < n; i++) { // iterate through all values in rt and out
      t = rt[i] - t0; // subtract non-decisison time from the response time
      if (t <= 0) {
        out[i] = 0; // if t=0, return 0 instead of NaN
        continue;
      }

      ks = ks_Nav(t / (a * a), eps); // number of terms in small time sum approximation
      kl = kl_Nav(t / (a * a), eps); // number of terms in sum approximation
      double log_t3 = 1.5 * log(t);
      if (response[i] == 1) { // if response is "upper" use alternate parameters
        if (ks < kl) { // small time needs fewer terms than lorge time
          mult = log_a - LOG_2PI_2 - log_t3 - vprime * a * wprime
               - vprime * vprime * t / 2;
          out[i] = mult + log(small_sum_2014(t, a, wprime, ks));
        } else { // large time needs fewer terms
          mult = LOG_PI - log_a2 - vprime * a * wprime - vprime * vprime * t / 2;
          out[i] = mult + log(large_sum_Nav(t, a, wprime, kl));
        }
      } else { // else response is "lower"
        if (ks < kl) { // small time needs fewer terms than lorge time
          mult = log_a - LOG_2PI_2 - log_t3 -v * a * w - v * v * t / 2;
          out[i] = mult + log(small_sum_2014(t, a, w, ks));
        } else { // large time needs fewer terms
          mult =  LOG_PI - log_a2 - v * a * w - v * v * t / 2;
          out[i] = mult + log(large_sum_Nav(t, a, w, kl));
        }
      }
    }
  }
  else { // use variable drift rate (changes mult) and only use small time
    double log_a = log(a);
    int ks;
    for (int i = 0; i < n; i++) { // iterate through all values in rt and out
      t = rt[i] - t0; // subtract non-decisison time from the response time
      if (t <= 0) {
        out[i] = 0; // if t=0, return 0 instead of NaN
        continue;
      }

      ks = ks_Nav(t / (a * a), eps); // number of terms in sum approximation
      double log_t3 = 1.5 * log(t);
      double log_svt = 0.5 * log(1 + sv * t);
      if (response[i] == 1) { // if response is "upper" use alternate parameters
        mult = log_a - log_t3 - log_svt + (sv * a * a * wprime * wprime
             - 2 * vprime * a * wprime - vprime * vprime * t) / (2 + 2 * sv * t);
        out[i] = mult + log(small_sum_2014(t, a, wprime, ks));
      } else { // else response is "lower"
        mult = log_a - log_t3 - log_svt + (sv * a * a * w * w
             - 2 * v * a * w - v * v * t) / (2 + 2 * sv * t);
        out[i] = mult + log(small_sum_2014(t, a, w, ks));
      }
    }
  }
  return out;
}


// ks = BGK2014, 2017 style sum approximation, kl = Navarro2009
Rcpp::NumericVector fb_BGK_Nav_2017_log(const Rcpp::NumericVector& rt,
                                        Rcpp::LogicalVector response,
                                        const double& a, const double& v,
                                        const double& t0, const double& w,
                                        const double& sv, const double& eps)
{
  int n = rt.length(); // get number of response times

  if (response.length() != n) { // create valid binary response vector
    bool first_resp = response[0];
    Rcpp::LogicalVector resp(n, first_resp);
    response = resp;
  }

  Rcpp::NumericVector out(n);
  double t, mult;
  double vprime = -v;
  double wprime = 1 - w;

  if (sv < SV_THRESH) { // set sv=0 and use constant drift rate methods
    double log_a = log(a);
    double log_a2 = 2 * log_a;
    int ks, kl;
    for (int i = 0; i < n; i++) { // iterate through all values in rt and out
      t = rt[i] - t0; // subtract non-decisison time from the response time
      if (t <= 0) {
        out[i] = 0; // if t=0, return 0 instead of NaN
        continue;
      }

      kl = kl_Nav(t / (a * a), eps); // number of terms in sum approximation
      double log_t3 = 1.5 * log(t);
      if (response[i] == 1) { // if response is "upper" use alternate parameters
        ks = ks_BGK(t / (a * a), wprime, eps); // number of terms in small time appx
        if (ks < kl) { // small time needs fewer terms than lorge time
          mult = log_a - LOG_2PI_2 - log_t3 - vprime * a * wprime
               - vprime * vprime * t / 2;
          out[i] = mult + log(small_sum_2017(t, a, wprime, ks));
        } else { // large time needs fewer terms
          mult = LOG_PI - log_a2 - vprime * a * wprime - vprime * vprime * t / 2;
          out[i] = mult + log(large_sum_Nav(t, a, wprime, kl));
        }
      } else { // else response is "lower"
        ks = ks_BGK(t / (a * a), w, eps); // number of terms in small time appx
        if (ks < kl) { // small time needs fewer terms than lorge time
          mult = log_a - LOG_2PI_2 - log_t3 -v * a * w - v * v * t / 2;
          out[i] = mult + log(small_sum_2017(t, a, w, ks));
        } else { // large time needs fewer terms
          mult =  LOG_PI - log_a2 - v * a * w - v * v * t / 2;
          out[i] = mult + log(large_sum_Nav(t, a, w, kl));
        }
      }
    }
  }
  else { // use variable drift rate (changes mult) and only use small time
    double log_a = log(a);
    int ks;
    for (int i = 0; i < n; i++) { // iterate through all values in rt and out
      t = rt[i] - t0; // subtract non-decisison time from the response time
      if (t <= 0) {
        out[i] = 0; // if t=0, return 0 instead of NaN
        continue;
      }

      double log_t3 = 1.5 * log(t);
      double log_svt = 0.5 * log(1 + sv * t);
      if (response[i] == 1) { // if response is "upper" use alternate parameters
        ks = ks_BGK(t / (a * a), wprime, eps); // number of terms in small time appx
        mult = log_a - log_t3 - log_svt + (sv * a * a * wprime * wprime
             - 2 * vprime * a * wprime - vprime * vprime * t) / (2 + 2 * sv * t);
        out[i] = mult + log(small_sum_2017(t, a, wprime, ks));
      } else { // else response is "lower"
        ks = ks_BGK(t / (a * a), w, eps); // number of terms in small time appx
        mult = log_a - log_t3 - log_svt + (sv * a * a * w * w
             - 2 * v * a * w - v * v * t) / (2 + 2 * sv * t);
        out[i] = mult + log(small_sum_2017(t, a, w, ks));
      }
    }
  }
  return out;
}


// ks = BGK2014, 2014 style sum approximation, kl = Navarro2009
Rcpp::NumericVector fb_BGK_Nav_2014_log(const Rcpp::NumericVector& rt,
                                        Rcpp::LogicalVector response,
                                        const double& a, const double& v,
                                        const double& t0, const double& w,
                                        const double& sv, const double& eps)
{
  int n = rt.length(); // get number of response times

  if (response.length() != n) { // create valid binary response vector
    bool first_resp = response[0];
    Rcpp::LogicalVector resp(n, first_resp);
    response = resp;
  }

  Rcpp::NumericVector out(n);
  double t, mult;
  double vprime = -v;
  double wprime = 1 - w;

  if (sv < SV_THRESH) { // set sv=0 and use constant drift rate methods
    double log_a = log(a);
    double log_a2 = 2 * log_a;
    int ks, kl;
    for (int i = 0; i < n; i++) { // iterate through all values in rt and out
      t = rt[i] - t0; // subtract non-decisison time from the response time
      if (t <= 0) {
        out[i] = 0; // if t=0, return 0 instead of NaN
        continue;
      }

      kl = kl_Nav(t / (a * a), eps); // number of terms in sum approximation
      double log_t3 = 1.5 * log(t);
      if (response[i] == 1) { // if response is "upper" use alternate parameters
        ks = ks_BGK(t / (a * a), wprime, eps); // number of terms in small time appx
        if (ks < kl) { // small time needs fewer terms than lorge time
          mult = log_a - LOG_2PI_2 - log_t3 - vprime * a * wprime
               - vprime * vprime * t / 2;
          out[i] = mult + log(small_sum_2014(t, a, wprime, ks));
        } else { // large time needs fewer terms
          mult = LOG_PI - log_a2 - vprime * a * wprime - vprime * vprime * t / 2;
          out[i] = mult + log(large_sum_Nav(t, a, wprime, kl));
        }
      } else { // else response is "lower"
        ks = ks_BGK(t / (a * a), w, eps); // number of terms in small time appx
        if (ks < kl) { // small time needs fewer terms than lorge time
          mult = log_a - LOG_2PI_2 - log_t3 -v * a * w - v * v * t / 2;
          out[i] = mult + log(small_sum_2014(t, a, w, ks));
        } else { // large time needs fewer terms
          mult =  LOG_PI - log_a2 - v * a * w - v * v * t / 2;
          out[i] = mult + log(large_sum_Nav(t, a, w, kl));
        }
      }
    }
  }
  else { // use variable drift rate (changes mult) and only use small time
    double log_a = log(a);
    int ks;
    for (int i = 0; i < n; i++) { // iterate through all values in rt and out
      t = rt[i] - t0; // subtract non-decisison time from the response time
      if (t <= 0) {
        out[i] = 0; // if t=0, return 0 instead of NaN
        continue;
      }

      double log_t3 = 1.5 * log(t);
      double log_svt = 0.5 * log(1 + sv * t);
      if (response[i] == 1) { // if response is "