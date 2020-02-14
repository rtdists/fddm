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
    double log_a = log(a);
    double log_t3 = 1.5 * log(t);
    if (sv < SV_THRESH) { // no sv
      mult = log_a - LOG_2PI_2 - log_t3 - v * a * w - v * v * t / 2;
      return mult + log(summ(t, a, w, ks, eps / exp(mult)));
    } else { // sv
      double log_svt = 0.5 * log(1 + sv * t);
      mult = log_a - log_t3 - log_svt + (sv * a * a * w * w
            - 2 * v * a * w - v * v * t) / (2 + 2 * sv * t);
      return mult + log(summ(t, a, w, ks, eps / exp(mult)));
    }
  } else { // no log
    if (sv < SV_THRESH) { // no sv
      mult = a * exp(-v * a * w - v * v * t / 2) / sqrt(2 * M_PI * t * t * t);
      return mult * summ(t, a, w, ks, eps / mult);
    } else { // sv
      mult = a * exp((-v * v * t - 2 * v * a * w + sv * a * a * w * w)
            / (2 + 2 * sv * t)) / sqrt(t * t * t + sv * t * t * t * t);
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
  // note: sv is not used because there is only a constant drift rate variant
  // also numm and summ are not used because there is only one option for each
  int kl = kl_Nav(t / (a * a), w, eps);
  double mult;
  
  if (log_prob) { // log
    double log_a2 = 2 * log(a);
    mult =  LOG_PI - log_a2 - v * a * w - v * v * t / 2;
    return mult + log(large_sum_Nav(t, a, w, kl, eps));
  } else { // no log
    mult = M_PI * exp(-v * a * w - v * v * t / 2) / (a * a);
    return mult * large_sum_Nav(t, a, w, kl, eps);
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
      double log_a = log(a);
      double log_t3 = 1.5 * log(t);
      if (sv < SV_THRESH) { // no sv
        mult = log_a - LOG_2PI_2 - log_t3 - v * a * w - v * v * t / 2;
        return mult + log(summ(t, a, w, ks, eps));
      } else { // sv
        double log_svt = 0.5 * log(1 + sv * t);
        mult = log_a - log_t3 - log_svt + (sv * a * a * w * w
              - 2 * v * a * w - v * v * t) / (2 + 2 * sv * t);
        return mult + log(summ(t, a, w, ks, eps));
      }
    } else { // no log
      if (sv < SV_THRESH) { // no sv
        mult = a * exp(-v * a * w - v * v * t / 2) / sqrt(2 * M_PI * t * t * t);
        return mult * summ(t, a, w, ks, eps);
      } else { // sv
        mult = a * exp((-v * v * t - 2 * v * a * w + sv * a * a * w * w)
              / (2 + 2 * sv * t)) / sqrt(t * t * t + sv * t * t * t * t);
        return mult * summ(t, a, w, ks, eps);
      }
    }
  } else { // large time appx is more efficient
    if (log_prob) { // log
      double log_a2 = 2 * log(a);
      mult =  LOG_PI - log_a2 - v * a * w - v * v * t / 2;
      return mult + log(large_sum_Nav(t, a, w, kl, eps));
    } else { // no log
      mult = M_PI * exp(-v * a * w - v * v * t / 2) / (a * a);
      return mult * large_sum_Nav(t, a, w, kl, eps);
    }
  }
}
