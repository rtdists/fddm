// Functions to evaluate the DDM PDF for specific criteria

#include "funcs.h"




//////////                                                            //////////
///////////////////////////////// Small Time ///////////////////////////////////
//////////                                                            //////////

double fs(const double& t, const double& a, const double& v, const double& w,
          const double& sv, const bool& log_prob, const double& eps,
          NummFunc numm, SummFunc summ)
{
  int ks = 0;
  if (numm != nullptr) { // don't calculate ks for Foster-style
    ks = numm(t / (a * a), w, eps);
  }
  double mult;

  if (log_prob) { // log
    if (sv < SV_THRESH) { // no sv
      mult = log(a) - LOG_2PI_2 - 1.5 * log(t) - v * a * w - v * v * t / 2;
      return mult + log(summ(t, a, w, ks, eps / exp(mult)));
    } else { // sv
      mult = log(a) - 1.5 * log(t) - 0.5 * log(1 + sv * t)
            + (sv * a * a * w * w - 2 * v * a * w - v * v * t)
            / (2 + 2 * sv * t);
      return mult + log(summ(t, a, w, ks, eps / exp(mult)));
    }
  } else { // no log
    if (sv < SV_THRESH) { // no sv
      mult = a * exp(-v * a * w - v * v * t / 2) / (t * sqrt(2 * M_PI * t));
      return mult * summ(t, a, w, ks, eps / mult);
    } else { // sv
      mult = a * exp((sv * a * a * w * w - 2 * v * a * w - v * v * t)
            / (2 + 2 * sv * t)) / (t * sqrt(t + sv * t * t));
      return mult * summ(t, a, w, ks, eps / mult);
    }
  }
}



//////////                                                            //////////
///////////////////////////////// Large Time ///////////////////////////////////
//////////                                                            //////////

double fl(const double& t, const double& a, const double& v, const double& w,
          const double& sv, const bool& log_prob, const double& eps,
          NummFunc numm, SummFunc summ)
{
  // also numm and summ are not used because there is only one option for each
  int kl = kl_Nav(t / (a * a), w, eps);
  double mult;
  
  if (log_prob) { // log
    if (sv < SV_THRESH) { // no sv
      mult =  LOG_PI - 2 * log(a) - v * a * w - v * v * t / 2;
      return mult + log(large_sum_Nav(t, a, w, kl, eps));
    } else { // sv
      mult = LOG_PI + LOG_2PI_2 - 2 * log(a) - 0.5 * log(1 + sv * t)
            + (sv * a * a * w * w - 2 * v * a * w - v * v * t)
            / (2 + 2 * sv * t);
      return mult + log(large_sum_Nav(t, a, w, kl, eps));
    }
  } else { // no log
    if (sv < SV_THRESH) { // no sv
      mult = M_PI * exp(-v * a * w - v * v * t / 2) / (a * a);
      return mult * large_sum_Nav(t, a, w, kl, eps);
    } else { // sv
      mult = M_PI * sqrt(2 * M_PI / (1 + sv * t))
            * exp((sv * a * a * w * w - 2 * v * a * w - v * v * t)
                  / (2 + 2 * sv * t)) / (a * a);
      return mult * large_sum_Nav(t, a, w, kl, eps);
    }
  }
}



//////////                                                            //////////
///////////////////////////////// Both Times ///////////////////////////////////
//////////                                                            //////////

double fb(const double& t, const double& a, const double& v, const double& w,
          const double& sv, const bool& log_prob, const double& eps,
          NummFunc numm, SummFunc summ)
{
  int ks = numm(t / (a * a), w, eps);
  int kl = kl_Nav(t / (a * a), w, eps);
  double mult;
  
  if (ks < kl) { // small time appx is more efficient
    if (log_prob) { // log
      if (sv < SV_THRESH) { // no sv
        mult = log(a) - LOG_2PI_2 - 1.5 * log(t) - v * a * w - v * v * t / 2;
        return mult + log(summ(t, a, w, ks, eps / exp(mult)));
      } else { // sv
        mult = log(a) - 1.5 * log(t) - 0.5 * log(1 + sv * t)
        + (sv * a * a * w * w - 2 * v * a * w - v * v * t)
        / (2 + 2 * sv * t);
        return mult + log(summ(t, a, w, ks, eps / exp(mult)));
      }
    } else { // no log
      if (sv < SV_THRESH) { // no sv
        mult = a * exp(-v * a * w - v * v * t / 2) / (t * sqrt(2 * M_PI * t));
        return mult * summ(t, a, w, ks, eps);
      } else { // sv
        mult = a * exp((sv * a * a * w * w - 2 * v * a * w - v * v * t)
              / (2 + 2 * sv * t)) / (t * sqrt(t + sv * t * t));
        return mult * summ(t, a, w, ks, eps);
      }
    }
  } else { // large time appx is more efficient
    if (log_prob) { // log
      if (sv < SV_THRESH) { // no sv
        mult =  LOG_PI - 2 * log(a) - v * a * w - v * v * t / 2;
        return mult + log(large_sum_Nav(t, a, w, kl, eps));
      } else { // sv
        mult = LOG_PI + LOG_2PI_2 - 2 * log(a) - 0.5 * log(1 + sv * t)
        + (sv * a * a * w * w - 2 * v * a * w - v * v * t)
        / (2 + 2 * sv * t);
        return mult + log(large_sum_Nav(t, a, w, kl, eps));
      }
    } else { // no log
      if (sv < SV_THRESH) { // no sv
        mult = M_PI * exp(-v * a * w - v * v * t / 2) / (a * a);
        return mult * large_sum_Nav(t, a, w, kl, eps);
      } else { // sv
        mult = M_PI * sqrt(2 * M_PI / (1 + sv * t))
        * exp((sv * a * a * w * w - 2 * v * a * w - v * v * t)
                / (2 + 2 * sv * t)) / (a * a);
                return mult * large_sum_Nav(t, a, w, kl, eps);
      }
    }
  }
}
