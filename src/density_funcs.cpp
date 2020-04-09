// Functions to evaluate the DDM PDF for specific criteria

#include "funcs.hpp"




//////////                                                            //////////
///////////////////////////////// Small Time ///////////////////////////////////
//////////                                                            //////////

double fs(const double& t, const double& a, const double& v, const double& w,
          const double& sv, const bool& log_prob, const double& eps,
          NummFunc numm, SummFunc summ, int ks)
{
  // don't calculate ks for Foster-style
  if (numm != nullptr && ks == -1) { // ks = -1: default argument from cpp_dfddm
    ks = numm(t / (a*a), w, eps);
  }
  double mult;

  if (log_prob) { // log
    if (sv < SV_THRESH) { // no sv
      mult = log(a) - LOG_2PI_2 - 1.5 * log(t) - v * a * w - v*v * t / 2;
      return mult + log(summ(t, a, w, ks, eps / exp(mult)));
    } else { // sv
      mult = log(a) - 1.5 * log(t) - LOG_2PI_2 - 0.5 * log(1 + sv*sv * t)
            + (sv*sv * a*a * w*w - 2 * v * a * w - v*v * t)
            / (2 + 2 * sv*sv * t);
      return mult + log(summ(t, a, w, ks, eps / exp(mult)));
    }
  } else { // no log
    if (sv < SV_THRESH) { // no sv
      mult = a * exp(-v * a * w - v * v * t / 2) / (t * SQRT_2PI * sqrt(t));
      return mult * summ(t, a, w, ks, eps / mult);
    } else { // sv
      mult = a * exp((sv*sv * a*a * w*w - 2 * v * a * w - v*v * t)
            / (2 + 2 * sv*sv * t)) / (t * SQRT_2PI * sqrt(t + sv*sv * t*t));
      return mult * summ(t, a, w, ks, eps / mult);
    }
  }
}



//////////                                                            //////////
///////////////////////////////// Large Time ///////////////////////////////////
//////////                                                            //////////

double fl(const double& t, const double& a, const double& v, const double& w,
          const double& sv, const bool& log_prob, const double& eps,
          NummFunc numm, SummFunc summ, int kl)
{
  // note: numm and summ are not used because there is only one option for each
  if (kl == -1) { // kl = -1: default argument from cpp_dfddm
    kl = kl_Nav(t / (a*a), w, eps);
  }
  double mult;

  if (log_prob) { // log
    if (sv < SV_THRESH) { // no sv
      mult =  LOG_PI - 2 * log(a) - v * a * w - v * v * t / 2;
      return mult + log(large_sum_Nav(t, a, w, kl, eps));
    } else { // sv
      mult = LOG_PI - 2 * log(a) - 0.5 * log(1 + sv*sv * t)
            + (sv*sv * a*a * w*w - 2 * v * a * w - v*v * t)
            / (2 + 2 * sv*sv * t);
      return mult + log(large_sum_Nav(t, a, w, kl, eps));
    }
  } else { // no log
    if (sv < SV_THRESH) { // no sv
      mult = M_PI * exp(-v * a * w - v*v * t / 2) / (a*a);
      return mult * large_sum_Nav(t, a, w, kl, eps);
    } else { // sv
      mult = M_PI * exp((sv*sv * a*a * w*w - 2 * v * a * w - v*v * t)
            / (2 + 2 * sv*sv * t)) / (a*a * sqrt(1 + sv*sv * t));
      return mult * large_sum_Nav(t, a, w, kl, eps);
    }
  }
}



//////////                                                            //////////
///////////////////////////////// Both Times ///////////////////////////////////
//////////                                                            //////////

double fb(const double& t, const double& a, const double& v, const double& w,
          const double& sv, const bool& log_prob, const double& eps,
          NummFunc numm, SummFunc summ, int k)
{
  // note: k isn't used
  int ks = numm(t / (a*a), w, eps);
  int kl = kl_Nav(t / (a*a), w, eps);

  if (ks < kl) { // small time appx is more efficient
    return fs(t, a, v, w, sv, log_prob, eps, numm, summ, ks);
  } else { // large time appx is more efficient
    return fl(t, a, v, w, sv, log_prob, eps, numm, summ, kl);
  }
}
