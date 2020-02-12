// Functions to evaluate the DDM PDF for specific criteria

#include <Rcpp.h>
// [[Rcpp::plugins(cpp11)]]
#define _USE_MATH_DEFINES
#include <math.h>
#include "funcs.h"




//////////                                                           //////////
///////////////////////////////// Small Time //////////////////////////////////
//////////                                                           //////////

// Use term < eps BGK2017 style sum approximation, with minimum terms
Rcpp::NumericVector fs_eps_2017(const Rcpp::NumericVector& rt,
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

  if (sv < SV_THRESH) { // set sv=0 and use constant drift rate method
    for (int i = 0; i < n; i++) { // iterate through all values in rt and out
      t = rt[i] - t0; // subtract non-decisison time from the response time
      if (t <= 0) {
        out[i] = 0; // if t=0, return 0 instead of NaN
        continue;
      }

      if (response[i] == 1) { // if response is "upper" use alternate parameters
        mult = a * exp(-vprime * a * wprime - vprime * vprime * t / 2)
             / sqrt(2 * M_PI * t * t * t);
        out[i] = mult * small_sum_eps_17(t, a, wprime, eps/mult);
      } else { // else response is "lower"
        mult = a * exp(-v * a * w - v * v * t / 2) / sqrt(2 * M_PI * t * t * t);
        out[i] = mult * small_sum_eps_17(t, a, w, eps/mult);
      }
    }
  } else { // use variable drift rate (changes mult)
    for (int i = 0; i < n; i++) { // iterate through all values in rt and out
      t = rt[i] - t0; // subtract non-decisison time from the response time
      if (t <= 0) {
        out[i] = 0; // if t=0, return 0 instead of NaN
        continue;
      }

      if (response[i] == 1) { // if response is "upper" use alternate parameters
        mult = a * exp((-vprime * vprime * t - 2 * vprime * a * wprime
                        + sv * a * a * wprime * wprime) / (2 + 2 * sv * t))
             / sqrt(t * t * t + sv * t * t * t * t);
        out[i] = mult * small_sum_eps_17(t, a, wprime, eps/mult);
      } else { // else response is "lower"
        mult = a * exp((-v * v * t - 2 * v * a * w
                        + sv * a * a * w * w) / (2 + 2 * sv * t))
             / sqrt(t * t * t + sv * t * t * t * t);
        out[i] = mult * small_sum_eps_17(t, a, w, eps/mult);
      }
    }
  }
  return out;
}


// Use term < eps BGK2014 style sum approximation, with minimum terms
Rcpp::NumericVector fs_eps_2014(const Rcpp::NumericVector& rt,
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

  if (sv < SV_THRESH) { // set sv=0 and use constant drift rate method
    for (int i = 0; i < n; i++) { // iterate through all values in rt and out
      t = rt[i] - t0; // subtract non-decisison time from the response time
      if (t <= 0) {
        out[i] = 0; // if t=0, return 0 instead of NaN
        continue;
      }

      if (response[i] == 1) { // if response is "upper" use alternate parameters
        mult = a * exp(-vprime * a * wprime - vprime * vprime * t / 2)
             / sqrt(2 * M_PI * t * t * t);
        out[i] = mult * small_sum_eps_14(t, a, wprime, eps/mult);
      } else { // else response is "lower"
        mult = a * exp(-v * a * w - v * v * t / 2) / sqrt(2 * M_PI * t * t * t);
        out[i] = mult * small_sum_eps_14(t, a, w, eps/mult);
      }
    }
  } else { // use variable drift rate (changes mult)
    for (int i = 0; i < n; i++) { // iterate through all values in rt and out
      t = rt[i] - t0; // subtract non-decisison time from the response time
      if (t <= 0) {
        out[i] = 0; // if t=0, return 0 instead of NaN
        continue;
      }

      if (response[i] == 1) { // if response is "upper" use alternate parameters
        mult = a * exp((-vprime * vprime * t - 2 * vprime * a * wprime
                        + sv * a * a * wprime * wprime) / (2 + 2 * sv * t))
             / sqrt(t * t * t + sv * t * t * t * t);
        out[i] = mult * small_sum_eps_14(t, a, wprime, eps/mult);
      } else { // else response is "lower"
        mult = a * exp((-v * v * t - 2 * v * a * w
                        + sv * a * a * w * w) / (2 + 2 * sv * t))
             / sqrt(t * t * t + sv * t * t * t * t);
        out[i] = mult * small_sum_eps_14(t, a, w, eps/mult);
      }
    }
  }
  return out;
}


// Use Navarro2009 number of terms for 2017 style sum approximation
Rcpp::NumericVector fs_Nav_2017(const Rcpp::NumericVector& rt,
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

  if (sv < SV_THRESH) { // set sv=0 and use constant drift rate method
    for (int i = 0; i < n; i++) { // iterate through all values in rt and out
      t = rt[i] - t0; // subtract non-decisison time from the response time
      if (t <= 0) {
        out[i] = 0; // if t=0, return 0 instead of NaN
        continue;
      }

      ks = ks_Nav(t / (a * a), eps); // number of terms in sum approximation
      if (response[i] == 1) { // if response is "upper" use alternate parameters
        mult = a * exp(-vprime * a * wprime - vprime * vprime * t / 2)
             / sqrt(2 * M_PI * t * t * t);
        out[i] = mult * small_sum_2017(t, a, wprime, ks);
      } else { // else response is "lower"
        mult = a * exp(-v * a * w - v * v * t / 2) / sqrt(2 * M_PI * t * t * t);
        out[i] = mult * small_sum_2017(t, a, w, ks);
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
      if (response[i] == 1) { // if response is "upper" use alternate parameters
        mult = a * exp((-vprime * vprime * t - 2 * vprime * a * wprime
                        + sv * a * a * wprime * wprime) / (2 + 2 * sv * t))
             / sqrt(t * t * t + sv * t * t * t * t);
        out[i] = mult * small_sum_2017(t, a, wprime, ks);
      } else { // else response is "lower"
        mult = a * exp((-v * v * t - 2 * v * a * w
                        + sv * a * a * w * w) / (2 + 2 * sv * t))
             / sqrt(t * t * t + sv * t * t * t * t);
        out[i] = mult * small_sum_2017(t, a, w, ks);
      }
    }
  }
  return out;
}


// Use Navarro2009 number of terms for 2014 style sum approximation
Rcpp::NumericVector fs_Nav_2014(const Rcpp::NumericVector& rt,
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

  if (sv < SV_THRESH) { // set sv=0 and use constant drift rate method
    for (int i = 0; i < n; i++) { // iterate through all values in rt and out
      t = rt[i] - t0; // subtract non-decisison time from the response time
      if (t <= 0) {
        out[i] = 0; // if t=0, return 0 instead of NaN
        continue;
      }

      ks = ks_Nav(t / (a * a), eps); // number of terms in sum approximation
      if (response[i] == 1) { // if response is "upper" use alternate parameters
        mult = a * exp(-vprime * a * wprime - vprime * vprime * t / 2)
             / sqrt(2 * M_PI * t * t * t);
        out[i] = mult * small_sum_2014(t, a, wprime, ks);
      } else { // else response is "lower"
        mult = a * exp(-v * a * w - v * v * t / 2) / sqrt(2 * M_PI * t * t * t);
        out[i] = mult * small_sum_2014(t, a, w, ks);
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
      if (response[i] == 1) { // if response is "upper" use alternate parameters
        mult = a * exp((-vprime * vprime * t - 2 * vprime * a * wprime
                        + sv * a * a * wprime * wprime) / (2 + 2 * sv * t))
             / sqrt(t * t * t + sv * t * t * t * t);
        out[i] = mult * small_sum_2014(t, a, wprime, ks);
      } else { // else response is "lower"
        mult = a * exp((-v * v * t - 2 * v * a * w
                        + sv * a * a * w * w) / (2 + 2 * sv * t))
             / sqrt(t * t * t + sv * t * t * t * t);
        out[i] = mult * small_sum_2014(t, a, w, ks);
      }
    }
  }
  return out;
}


// Use BGK2014 number of terms for 2017 style sum approximation
Rcpp::NumericVector fs_BGK_2017(const Rcpp::NumericVector& rt,
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

  if (sv < SV_THRESH) { // set sv=0 and use constant drift rate method
    for (int i = 0; i < n; i++) { // iterate through all values in rt and out
      t = rt[i] - t0; // subtract non-decisison time from the response time
      if (t <= 0) {
        out[i] = 0; // if t=0, return 0 instead of NaN
        continue;
      }

      if (response[i] == 1) { // if response is "upper" use alternate parameters
        ks = ks_BGK(t / (a * a), wprime, eps); // number of terms in sum approximation
        mult = a * exp(-vprime * a * wprime - vprime * vprime * t / 2)
             / sqrt(2 * M_PI * t * t * t);
        out[i] = mult * small_sum_2017(t, a, wprime, ks);
      } else { // else response is "lower"
        ks = ks_BGK(t / (a * a), w, eps); // number of terms in sum approximation
        mult = a * exp(-v * a * w - v * v * t / 2) / sqrt(2 * M_PI * t * t * t);
        out[i] = mult * small_sum_2017(t, a, w, ks);
      }
    }
  } else { // use variable drift rate (changes mult)
    for (int i = 0; i < n; i++) { // iterate through all values in rt and out
      t = rt[i] - t0; // subtract non-decisison time from the response time
      if (t <= 0) {
        out[i] = 0; // if t=0, return 0 instead of NaN
        continue;
      }

      if (response[i] == 1) { // if response is "upper" use alternate parameters
        ks = ks_BGK(t / (a * a), wprime, eps); // number of terms in sum approximation
        mult = a * exp((-vprime * vprime * t - 2 * vprime * a * wprime
                        + sv * a * a * wprime * wprime) / (2 + 2 * sv * t))
             / sqrt(t * t * t + sv * t * t * t * t);
        out[i] = mult * small_sum_2017(t, a, wprime, ks);
      } else { // else response is "lower"
        ks = ks_BGK(t / (a * a), w, eps); // number of terms in sum approximation
        mult = a * exp((-v * v * t - 2 * v * a * w
                        + sv * a * a * w * w) / (2 + 2 * sv * t))
             / sqrt(t * t * t + sv * t * t * t * t);
        out[i] = mult * small_sum_2017(t, a, w, ks);
      }
    }
  }
  return out;
}


// Use BGK2014 number of terms for 2014 style sum approximation
Rcpp::NumericVector fs_BGK_2014(const Rcpp::NumericVector& rt,
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

  if (sv < SV_THRESH) { // set sv=0 and use constant drift rate method
    for (int i = 0; i < n; i++) { // iterate through all values in rt and out
      t = rt[i] - t0; // subtract non-decisison time from the response time
      if (t <= 0) {
        out[i] = 0; // if t=0, return 0 instead of NaN
        continue;
      }

      if (response[i] == 1) { // if response is "upper" use alternate parameters
        ks = ks_BGK(t / (a * a), wprime, eps); // number of terms in sum approximation
        mult = a * exp(-vprime * a * wprime - vprime * vprime * t / 2)
             / sqrt(2 * M_PI * t * t * t);
        out[i] = mult * small_sum_2014(t, a, wprime, ks);
      } else { // else response is "lower"
        ks = ks_BGK(t / (a * a), w, eps); // number of terms in sum approximation
        mult = a * exp(-v * a * w - v * v * t / 2) / sqrt(2 * M_PI * t * t * t);
        out[i] = mult * small_sum_2014(t, a, w, ks);
      }
    }
  } else { // use variable drift rate (changes mult)
    for (int i = 0; i < n; i++) { // iterate through all values in rt and out
      t = rt[i] - t0; // subtract non-decisison time from the response time
      if (t <= 0) {
        out[i] = 0; // if t=0, return 0 instead of NaN
        continue;
      }

      if (response[i] == 1) { // if response is "upper" use alternate parameters
        ks = ks_BGK(t / (a * a), wprime, eps); // number of terms in sum approximation
        mult = a * exp((-vprime * vprime * t - 2 * vprime * a * wprime
                        + sv * a * a * wprime * wprime) / (2 + 2 * sv * t))
             / sqrt(t * t * t + sv * t * t * t * t);
        out[i] = mult * small_sum_2014(t, a, wprime, ks);
      } else { // else response is "lower"
        ks = ks_BGK(t / (a * a), w, eps); // number of terms in sum approximation
        mult = a * exp((-v * v * t - 2 * v * a * w
                        + sv * a * a * w * w) / (2 + 2 * sv * t))
             / sqrt(t * t * t + sv * t * t * t * t);
        out[i] = mult * small_sum_2014(t, a, w, ks);
      }
    }
  }
  return out;
}



//////////                                                           //////////
///////////////////////////////// Large Time //////////////////////////////////
//////////                                                           //////////

// Use Navarro2009 number of terms for sum approximation
Rcpp::NumericVector fl_Nav_2009(const Rcpp::NumericVector& rt,
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

  for (int i = 0; i < n; i++) { // iterate through all values in rt and out
    t = rt[i] - t0; // subtract non-decisison time from the response time
    if (t <= 0) {
      out[i] = 0; // if t=0, return 0 instead of NaN
      continue;
    }

    kl = kl_Nav(t / (a * a), eps); // number of terms in sum approximation
    if (response[i] == 1) { // if response is "upper" use alternate parameters
      mult = M_PI * exp(-vprime * a * wprime - vprime * vprime * t / 2)
           / (a * a);
      out[i] = mult * large_sum_Nav(t, a, wprime, kl);
    } else { // else response is "lower"
      mult = M_PI * exp(-v * a * w - v * v * t / 2) / (a * a);
      out[i] = mult * large_sum_Nav(t, a, w, kl);
    }
  }
  return out;
}



//////////                                                           //////////
///////////////////////////////// Both Times //////////////////////////////////
//////////                                                           //////////

// ks = Navarro2009, 2017 style sum approximation, kl = Navarro2009
Rcpp::NumericVector fb_Nav_Nav_2017(const Rcpp::NumericVector& rt,
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
    int ks, kl;
    for (int i = 0; i < n; i++) { // iterate through all values in rt and out
      t = rt[i] - t0; // subtract non-decisison time from the response time
      if (t <= 0) {
        out[i] = 0; // if t=0, return 0 instead of NaN
        continue;
      }

      ks = ks_Nav(t / (a * a), eps); // number of terms in small time sum approximation
      kl = kl_Nav(t / (a * a), eps); // number of terms in sum approximation
      if (response[i] == 1) { // if response is "upper" use alternate parameters
        if (ks < kl) { // small time needs fewer terms than lorge time
          mult = a * exp((-vprime * vprime * t - 2 * vprime * a * wprime
                          + sv * a * a * wprime * wprime) / (2 + 2 * sv * t))
               / sqrt(t * t * t + sv * t * t * t * t);
          out[i] = mult * small_sum_2017(t, a, wprime, ks);
        } else { // large time needs fewer terms
          mult = M_PI * exp(-vprime * a * wprime - vprime * vprime * t / 2)
               / (a * a);
          out[i] = mult * large_sum_Nav(t, a, wprime, kl);
        }
      } else { // else response is "lower"
        if (ks < kl) { // small time needs fewer terms than lorge time
          mult = a * exp((-v * v * t - 2 * v * a * w
                          + sv * a * a * w * w) / (2 + 2 * sv * t))
               / sqrt(t * t * t + sv * t * t * t * t);
          out[i] = mult * small_sum_2017(t, a, w, ks);
        } else { // large time needs fewer terms
          mult = M_PI * exp(-v * a * w - v * v * t / 2) / (a * a);
          out[i] = mult * large_sum_Nav(t, a, w, kl);
        }
      }
    }
  } else { // use variable drift rate (changes mult) and only use small time
    int ks;
    for (int i = 0; i < n; i++) { // iterate through all values in rt and out
      t = rt[i] - t0; // subtract non-decisison time from the response time
      if (t <= 0) {
        out[i] = 0; // if t=0, return 0 instead of NaN
        continue;
      }

      ks = ks_Nav(t / (a * a), eps); // number of terms in sum approximation
      if (response[i] == 1) { // if response is "upper" use alternate parameters
        mult = a * exp(-vprime * a * wprime - vprime * vprime * t / 2)
             / sqrt(2 * M_PI * t * t * t);
        out[i] = mult * small_sum_2017(t, a, wprime, ks);
      } else { // else response is "lower"
        mult = a * exp(-v * a * w - v * v * t / 2) / sqrt(2 * M_PI * t * t * t);
        out[i] = mult * small_sum_2017(t, a, w, ks);
      }
    }
  }
  return out;
}


// ks = Navarro2009, 2014 style sum approximation, kl = Navarro2009
Rcpp::NumericVector fb_Nav_Nav_2014(const Rcpp::NumericVector& rt,
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
    int ks, kl;
    for (int i = 0; i < n; i++) { // iterate through all values in rt and out
      t = rt[i] - t0; // subtract non-decisison time from the response time
      if (t <= 0) {
        out[i] = 0; // if t=0, return 0 instead of NaN
        continue;
      }

      ks = ks_Nav(t / (a * a), eps); // number of terms in small time sum approximation
      kl = kl_Nav(t / (a * a), eps); // number of terms in sum approximation
      if (response[i] == 1) { // if response is "upper" use alternate parameters
        if (ks < kl) { // small time needs fewer terms than lorge time
          mult = a * exp((-vprime * vprime * t - 2 * vprime * a * wprime
                          + sv * a * a * wprime * wprime) / (2 + 2 * sv * t))
               / sqrt(t * t * t + sv * t * t * t * t);
          out[i] = mult * small_sum_2014(t, a, wprime, ks);
        } else { // large time needs fewer terms
          mult = M_PI * exp(-vprime * a * wprime - vprime * vprime * t / 2)
               / (a * a);
          out[i] = mult * large_sum_Nav(t, a, wprime, kl);
        }
      } else { // else response is "lower"
        if (ks < kl) { // small time needs fewer terms than lorge time
          mult = a * exp((-v * v * t - 2 * v * a * w
                          + sv * a * a * w * w) / (2 + 2 * sv * t))
               / sqrt(t * t * t + sv * t * t * t * t);
          out[i] = mult * small_sum_2014(t, a, w, ks);
        } else { // large time needs fewer terms
          mult = M_PI * exp(-v * a * w - v * v * t / 2) / (a * a);
          out[i] = mult * large_sum_Nav(t, a, w, kl);
        }
      }
    }
  } else { // use variable drift rate (changes mult) and only use small time
    int ks;
    for (int i = 0; i < n; i++) { // iterate through all values in rt and out
      t = rt[i] - t0; // subtract non-decisison time from the response time
      if (t <= 0) {
        out[i] = 0; // if t=0, return 0 instead of NaN
        continue;
      }

      ks = ks_Nav(t / (a * a), eps); // number of terms in sum approximation
      if (response[i] == 1) { // if response is "upper" use alternate parameters
        mult = a * exp(-vprime * a * wprime - vprime * vprime * t / 2)
             / sqrt(2 * M_PI * t * t * t);
        out[i] = mult * small_sum_2014(t, a, wprime, ks);
      } else { // else response is "lower"
        mult = a * exp(-v * a * w - v * v * t / 2) / sqrt(2 * M_PI * t * t * t);
        out[i] = mult * small_sum_2014(t, a, w, ks);
      }
    }
  }
  return out;
}


// ks = BGK2014, 2017 style sum approximation, kl = Navarro2009
Rcpp::NumericVector fb_BGK_Nav_2017(const Rcpp::NumericVector& rt,
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
    int ks, kl;
    for (int i = 0; i < n; i++) { // iterate through all values in rt and out
      t = rt[i] - t0; // subtract non-decisison time from the response time
      if (t <= 0) {
        out[i] = 0; // if t=0, return 0 instead of NaN
        continue;
      }

      kl = kl_Nav(t / (a * a), eps); // number of terms in sum approximation
      if (response[i] == 1) { // if response is "upper" use alternate parameters
        ks = ks_BGK(t / (a * a), wprime, eps); // number of terms in small time appx
        if (ks < kl) { // small time needs fewer terms than lorge time
          mult = a * exp((-vprime * vprime * t - 2 * vprime * a * wprime
                          + sv * a * a * wprime * wprime) / (2 + 2 * sv * t))
               / sqrt(t * t * t + sv * t * t * t * t);
          out[i] = mult * small_sum_2017(t, a, wprime, ks);
        } else { // large time needs fewer terms
          mult = M_PI * exp(-vprime * a * wprime - vprime * vprime * t / 2)
               / (a * a);
          out[i] = mult * large_sum_Nav(t, a, wprime, kl);
        }
      } else { // else response is "lower"
        ks = ks_BGK(t / (a * a), w, eps); // number of terms in small time appx
        if (ks < kl) { // small time needs fewer terms than lorge time
          mult = a * exp((-v * v * t - 2 * v * a * w
                          + sv * a * a * w * w) / (2 + 2 * sv * t))
               / sqrt(t * t * t + sv * t * t * t * t);
          out[i] = mult * small_sum_2017(t, a, w, ks);
        } else { // large time needs fewer terms
          mult = M_PI * exp(-v * a * w - v * v * t / 2) / (a * a);
          out[i] = mult * large_sum_Nav(t, a, w, kl);
        }
      }
    }
  } else { // use variable drift rate (changes mult) and only use small time
    int ks;
    for (int i = 0; i < n; i++) { // iterate through all values in rt and out
      t = rt[i] - t0; // subtract non-decisison time from the response time
      if (t <= 0) {
        out[i] = 0; // if t=0, return 0 instead of NaN
        continue;
      }

      if (response[i] == 1) { // if response is "upper" use alternate parameters
        ks = ks_BGK(t / (a * a), wprime, eps); // number of terms in small time appx
        mult = a * exp(-vprime * a * wprime - vprime * vprime * t / 2)
             / sqrt(2 * M_PI * t * t * t);
        out[i] = mult * small_sum_2017(t, a, wprime, ks);
      } else { // else response is "lower"
        ks = ks_BGK(t / (a * a), w, eps); // number of terms in small time appx
        mult = a * exp(-v * a * w - v * v * t / 2) / sqrt(2 * M_PI * t * t * t);
        out[i] = mult * small_sum_2017(t, a, w, ks);
      }
    }
  }
  return out;
}


// ks = BGK2014, 2014 style sum approximation, kl = Navarro2009
Rcpp::NumericVector fb_BGK_Nav_2014(const Rcpp::NumericVector& rt,
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
    int ks, kl;
    for (int i = 0; i < n; i++) { // iterate through all values in rt and out
      t = rt[i] - t0; // subtract non-decisison time from the response time
      if (t <= 0) {
        out[i] = 0; // if t=0, return 0 instead of NaN
        continue;
      }

      kl = kl_Nav(t / (a * a), eps); // number of terms in sum approximation
      if (response[i] == 1) { // if response is "upper" use alternate parameters
        ks = ks_BGK(t / (a * a), wprime, eps); // number of terms in small time appx
        if (ks < kl) { // small time needs fewer terms than lorge time
          mult = a * exp((-vprime * vprime * t - 2 * vprime * a * wprime
                          + sv * a * a * wprime * wprime) / (2 + 2 * sv * t))
               / sqrt(t * t * t + sv * t * t * t * t);
          out[i] = mult * small_sum_2014(t, a, wprime, ks);
        } else { // large time needs fewer terms
          mult = M_PI * exp(-vprime * a * wprime - vprime * vprime * t / 2)
               / (a * a);
          out[i] = mult * large_sum_Nav(t, a, wprime, kl);
        }
      } else { // else response is "lower"
        ks = ks_BGK(t / (a * a), w, eps); // number of terms in small time appx
        if (ks < kl) { // small time needs fewer terms than lorge time
          mult = a * exp((-v * v * t - 2 * v * a * w
                          + sv * a * a * w * w) / (2 + 2 * sv * t))
               / sqrt(t * t * t + sv * t * t * t * t);
          out[i] = mult * small_sum_2014(t, a, w, ks);
        } else { // large time needs fewer terms
          mult = M_PI * exp(-v * a * w - v * v * t / 2) / (a * a);
          out[i] = mult * large_sum_Nav(t, a, w, kl);
        }
      }
    }
  } else { // use variable drift rate (changes mult) and only use small time
    int ks;
    for (int i = 0; i < n; i++) { // iterate through all values in rt and out
      t = rt[i] - t0; // subtract non-decisison time from the response time
      if (t <= 0) {
        out[i] = 0; // if t=0, return 0 instead of NaN
        continue;
      }

      if (response[i] == 1) { // if response is "upper" use alternate parameters
        ks = ks_BGK(t / (a * a), wprime, eps); // number of terms in small time appx
        mult = a * exp(-vprime * a * wprime - vprime * vprime * t / 2)
             / sqrt(2 * M_PI * t * t * t);
        out[i] = mult * small_sum_2014(t, a, wprime, ks);
      } else { // else response is "lower"
        ks = ks_BGK(t / (a * a), w, eps); // number of terms in small time appx
        mult = a * exp(-v * a * w - v * v * t / 2) / sqrt(2 * M_PI * t * t * t);
        out[i] = mult * small_sum_2014(t, a, w, ks);
      }
    }
  }
  return out;
}
