// Functions to evaluate the DDM PDF for specific criteria

#include "funcs.h"




//////////                                                            //////////
///////////////////////////////// Small Time ///////////////////////////////////
//////////                                                            //////////

// Foster
double ff(const double& t, const double& a, const double& v,
          const double& w, const double& sv, const double& eps,
          NummFunc numm, SummFunc summ)
{
  // note: numm is not used
  double mult;

  if (sv < SV_THRESH) { // no sv
    mult = a * exp(-v * a * w - v * v * t / 2) / (t * SQRT_2PI * sqrt(t));
  } else { // sv
    mult = a * exp((sv*sv * a*a * w*w - 2 * v * a * w - v*v * t)
          / (2 + 2 * sv*sv * t)) / (t * SQRT_2PI * sqrt(t + sv*sv * t*t));
  }
  return mult * summ(t, a, w, 0, eps / mult);
}

double ff_log(const double& t, const double& a, const double& v,
              const double& w, const double& sv, const double& eps,
              NummFunc numm, SummFunc summ)
{
  // note: numm is not used
  double mult;

  if (sv < SV_THRESH) { // no sv
    mult = log(a) - LOG_2PI_2 - 1.5 * log(t) - v * a * w - v*v * t / 2;
  } else { // sv
    mult = log(a) - 1.5 * log(t) - LOG_2PI_2 - 0.5 * log(1 + sv*sv * t)
          + (sv*sv * a*a * w*w - 2 * v * a * w - v*v * t)
          / (2 + 2 * sv*sv * t);
  }
  return mult + log(summ(t, a, w, 0, eps / exp(mult)));
}

// Kesselmeier and Navarro
double fs(const double& t, const double& a, const double& v,
          const double& w, const double& sv, const double& eps,
          NummFunc numm, SummFunc summ)
{
  double mult_s;
  int ks;

  if (sv < SV_THRESH) { // no sv
    mult_s = exp(-v * a * w - v * v * t / 2);
  } else { // sv
    mult_s = exp((sv*sv * a*a * w*w - 2 * v * a * w - v*v * t)
                 / (2 + 2 * sv*sv * t)) / sqrt(1 + sv*sv * t);
  }
  ks = numm(t / (a*a), w, eps * a*a / mult_s);
  mult_s *= a / (t * SQRT_2PI * sqrt(t));
  return mult_s * summ(t, a, w, ks, 0.0);
}

double fs_log(const double& t, const double& a, const double& v,
              const double& w, const double& sv, const double& eps,
              NummFunc numm, SummFunc summ)
{
  double mult_s;
  int ks;

  if (sv < SV_THRESH) { // no sv
    mult_s = -v * a * w - v * v * t / 2;
  } else { // sv
    mult_s = (sv*sv * a*a * w*w - 2 * v * a * w - v*v * t)
             / (2 + 2 * sv*sv * t) - 0.5 * log(1 + sv*sv * t);
  }
  ks = numm(t / (a*a), w, eps * a*a / exp(mult_s));
  mult_s += log(a) - LOG_2PI_2 - 1.5 * log(t);
  return mult_s + log(summ(t, a, w, ks, 0.0));
}



//////////                                                            //////////
///////////////////////////////// Large Time ///////////////////////////////////
//////////                                                            //////////

double fl(const double& t, const double& a, const double& v,
          const double& w, const double& sv, const double& eps,
          NummFunc numm, SummFunc summ)
{
  // note: neither numm nor summ is used because there is only one option each
  double mult_l;
  int kl;

  if (sv < SV_THRESH) { // no sv
    mult_l = exp(-v * a * w - v*v * t / 2) / (a*a);
  } else { // sv
    mult_l = exp((sv*sv * a*a * w*w - 2 * v * a * w - v*v * t)
                 / (2 + 2 * sv*sv * t)) / (a*a * sqrt(1 + sv*sv * t));
  }
  kl = kl_Nav(t / (a*a), w, eps / mult_l);
  mult_l *= M_PI;
  return mult_l * large_sum_Nav(t, a, w, kl, 0.0);
}

double fl_log(const double& t, const double& a, const double& v,
              const double& w, const double& sv, const double& eps,
              NummFunc numm, SummFunc summ)
{
  // note: neither numm nor summ is used because there is only one option each
  double mult_l;
  int kl;

  if (sv < SV_THRESH) { // no sv
    mult_l = - v * a * w - v * v * t / 2 - 2 * log(a);
  } else { // sv
    mult_l = (sv*sv * a*a * w*w - 2 * v * a * w - v*v * t)
             / (2 + 2 * sv*sv * t) - 0.5 * log(1 + sv*sv * t) - 2 * log(a);
  }
  kl = kl_Nav(t / (a*a), w, eps / exp(mult_l));
  mult_l += LOG_PI;
  return mult_l + log(large_sum_Nav(t, a, w, kl, 0.0));
}



//////////                                                            //////////
///////////////////////////////// Both Times ///////////////////////////////////
//////////                                                            //////////

double fb(const double& t, const double& a, const double& v,
          const double& w, const double& sv, const double& eps,
          NummFunc numm, SummFunc summ)
{
  double mult_s, mult_l;
  int ks, kl;

  if (sv < SV_THRESH) { // no sv
    mult_s = exp(-v * a * w - v * v * t / 2);
    mult_l = exp(-v * a * w - v*v * t / 2) / (a*a);
  } else { // sv
    mult_s = exp((sv*sv * a*a * w*w - 2 * v * a * w - v*v * t)
                 / (2 + 2 * sv*sv * t)) / sqrt(1 + sv*sv * t);
    mult_l = exp((sv*sv * a*a * w*w - 2 * v * a * w - v*v * t)
                 / (2 + 2 * sv*sv * t)) / (a*a * sqrt(1 + sv*sv * t));
  }
  ks = numm(t / (a*a), w, eps * a*a / mult_s);
  mult_s *= a / (t * SQRT_2PI * sqrt(t));
  kl = kl_Nav(t / (a*a), w, eps / mult_l);
  mult_l *= M_PI;
  if (ks < kl) { // small-time is better
    return mult_s * summ(t, a, w, ks, 0.0);
  } else { // large-time is better
    return mult_l * large_sum_Nav(t, a, w, kl, 0.0);
  }
}

double fb_log(const double& t, const double& a, const double& v,
              const double& w, const double& sv, const double& eps,
              NummFunc numm, SummFunc summ)
{
  double mult_s, mult_l;
  int ks, kl;

  if (sv < SV_THRESH) { // no sv
    mult_s = -v * a * w - v * v * t / 2;
    mult_l = - v * a * w - v * v * t / 2 - 2 * log(a);
  } else { // sv
    mult_s = (sv*sv * a*a * w*w - 2 * v * a * w - v*v * t)
             / (2 + 2 * sv*sv * t) - 0.5 * log(1 + sv*sv * t);
    mult_l = (sv*sv * a*a * w*w - 2 * v * a * w - v*v * t)
             / (2 + 2 * sv*sv * t) - 0.5 * log(1 + sv*sv * t) - 2 * log(a);
  }
  ks = numm(t / (a*a), w, eps * a*a / exp(mult_s));
  kl = kl_Nav(t / (a*a), w, eps / exp(mult_l));
  mult_s += log(a) - LOG_2PI_2 - 1.5 * log(t);
  mult_l += LOG_PI;
  if (ks < kl) { // small-time is better
    return mult_s + log(summ(t, a, w, ks, 0.0));
  } else { // large-time is better
    return mult_l + log(large_sum_Nav(t, a, w, kl, 0.0));
  }
}
