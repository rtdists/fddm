// Functions to evaluate the DDM PDF for specific criteria



//////////                                                            //////////
///////////////////////////////// Small Time ///////////////////////////////////
//////////                                                            //////////


// Stop When Small Enough
double nff(const double& t, const double& a, const double& v,
          const double& w, const double& sv, const double& err,
          const double& sl_thresh, const NumFunc& numf, const SumFunc& sumf)
{ // note: sl_thresh, numf are not used
  double mult;

  if (sv <= SV_THRESH) { // no sv
    mult = a * exp(-v * a * w - v * v * t / 2) / (t * SQRT_2PI * sqrt(t));
  } else { // sv
    mult = a * exp((sv*sv * a*a * w*w - 2 * v * a * w - v*v * t)
                     / (2 + 2 * sv*sv * t)) / (t * SQRT_2PI * sqrt(t + sv*sv * t*t));
  }
  return mult * sumf(t, a, w, 0, err / mult);
}

double nff_log(const double& t, const double& a, const double& v,
              const double& w, const double& sv, const double& err,
              const double& sl_thresh,
              const NumFunc& numf, const SumFunc& sumf)
{ // note: sl_thresh, numf are not used
  double mult;

  if (sv <= SV_THRESH) { // no sv
    mult = log(a) - LOG_2PI_2 - 1.5 * log(t) - v * a * w - v*v * t / 2;
  } else { // sv
    mult = log(a) - 1.5 * log(t) - LOG_2PI_2 - 0.5 * log(1 + sv*sv * t)
    + (sv*sv * a*a * w*w - 2 * v * a * w - v*v * t)
    / (2 + 2 * sv*sv * t);
  }
  double temp = sumf(t, a, w, 0, err / exp(mult));
  if (temp > 0) {
    return mult + log(temp);
  } else{ // protect against -Inf
    return log(err) - LOG_100;
  }
}


// Gondan et al and Navarro et al
double nfs(const double& t, const double& a, const double& v,
          const double& w, const double& sv, const double& err,
          const double& sl_thresh, const NumFunc& numf, const SumFunc& sumf)
{ // note: sl_thresh is not used
  double mult_s;
  int ks;

  if (sv <= SV_THRESH) { // no sv
    mult_s = exp(-v * a * w - v * v * t / 2);
  } else { // sv
    mult_s = exp((sv*sv * a*a * w*w - 2 * v * a * w - v*v * t)
                   / (2 + 2 * sv*sv * t)) / sqrt(1 + sv*sv * t);
  }
  ks = numf(t / (a*a), w, err * a*a / mult_s);
  mult_s *= a / (t * SQRT_2PI * sqrt(t));
  return mult_s * sumf(t, a, w, ks, 0.0);
}

double nfs_log(const double& t, const double& a, const double& v,
              const double& w, const double& sv, const double& err,
              const double& sl_thresh,
              const NumFunc& numf, const SumFunc& sumf)
{ // note: sl_thresh is not used
  double mult_s;
  int ks;

  if (sv <= SV_THRESH) { // no sv
    mult_s = -v * a * w - v * v * t / 2;
  } else { // sv
    mult_s = (sv*sv * a*a * w*w - 2 * v * a * w - v*v * t)
    / (2 + 2 * sv*sv * t) - 0.5 * log(1 + sv*sv * t);
  }
  ks = numf(t / (a*a), w, err * a*a / exp(mult_s));
  mult_s += log(a) - LOG_2PI_2 - 1.5 * log(t);
  return mult_s + log(sumf(t, a, w, ks, 0.0));
}





//////////                                                            //////////
///////////////////////////////// Large Time ///////////////////////////////////
//////////                                                            //////////


// Navarro et al
double nfl(const double& t, const double& a, const double& v,
          const double& w, const double& sv, const double& err,
          const double& sl_thresh, const NumFunc& numf, const SumFunc& sumf)
{ // note: sl_thresh, numf, sumf are not used
  double mult_l;
  int kl;

  if (sv <= SV_THRESH) { // no sv
    mult_l = exp(-v * a * w - v*v * t / 2) / (a*a);
  } else { // sv
    mult_l = exp((sv*sv * a*a * w*w - 2 * v * a * w - v*v * t)
                   / (2 + 2 * sv*sv * t)) / (a*a * sqrt(1 + sv*sv * t));
  }
  kl = nkl_Nav(t / (a*a), w, err / mult_l);
  mult_l *= PI_CONST;
  return mult_l * nlarge_sum_Nav(t, a, w, kl, 0.0);
}

double nfl_log(const double& t, const double& a, const double& v,
              const double& w, const double& sv, const double& err,
              const double& sl_thresh,
              const NumFunc& numf, const SumFunc& sumf)
{ // note: sl_thresh, numf, sumf are not used
  double mult_l;
  int kl;

  if (sv <= SV_THRESH) { // no sv
    mult_l = - v * a * w - v * v * t / 2 - 2 * log(a);
  } else { // sv
    mult_l = (sv*sv * a*a * w*w - 2 * v * a * w - v*v * t)
    / (2 + 2 * sv*sv * t) - 0.5 * log(1 + sv*sv * t) - 2 * log(a);
  }
  kl = nkl_Nav(t / (a*a), w, err / exp(mult_l));
  mult_l += LOG_PI;
  return mult_l + log(nlarge_sum_Nav(t, a, w, kl, 0.0));
}



//////////                                                            //////////
///////////////////////////////// Both Times ///////////////////////////////////
//////////                                                            //////////


// Stop When Small Enough
double nfc(const double& t, const double& a, const double& v,
          const double& w, const double& sv, const double& err,
          const double& sl_thresh, const NumFunc& numf, const SumFunc& sumf)
{ // note: numf is not used
  double mult;

  if (t / (a*a) > sl_thresh) { // use large time
    if (sv <= SV_THRESH) { // no sv
      mult = exp(-v * a * w - v*v * t / 2) / (a*a);
    } else { // sv
      mult = exp((sv*sv * a*a * w*w - 2 * v * a * w - v*v * t)
                  / (2 + 2 * sv*sv * t)) / (a*a * sqrt(1 + sv*sv * t));
    }
    int kl = nkl_Nav(t / (a*a), w, err / mult);
    return PI_CONST * mult * nlarge_sum_Nav(t, a, w, kl, 0.0);
  } else { // use small time
    if (sv <= SV_THRESH) { // no sv
      mult = a * exp(-v * a * w - v * v * t / 2) / (t * SQRT_2PI * sqrt(t));
    } else { // sv
      mult = a * exp((sv*sv * a*a * w*w - 2 * v * a * w - v*v * t)
                       / (2 + 2 * sv*sv * t)) / (t * SQRT_2PI * sqrt(t + sv*sv * t*t));
    }
    return mult * sumf(t, a, w, 0, err / mult);
  }
}

double nfc_log(const double& t, const double& a, const double& v,
              const double& w, const double& sv, const double& err,
              const double& sl_thresh,
              const NumFunc& numf, const SumFunc& sumf)
{ // note: numf is not used
  double mult;

  if (t / (a*a) > sl_thresh) { // use large time
    if (sv <= SV_THRESH) { // no sv
      mult = - v * a * w - v * v * t / 2 - 2 * log(a);
    } else { // sv
      mult = (sv*sv * a*a * w*w - 2 * v * a * w - v*v * t)
      / (2 + 2 * sv*sv * t) - 0.5 * log(1 + sv*sv * t) - 2 * log(a);
    }
    int kl = nkl_Nav(t / (a*a), w, err / exp(mult));
    return LOG_PI + mult + log(nlarge_sum_Nav(t, a, w, kl, 0.0));
  } else { // use small time
    if (sv <= SV_THRESH) { // no sv
      mult = log(a) - LOG_2PI_2 - 1.5 * log(t) - v * a * w - v*v * t / 2;
    } else { // sv
      mult = log(a) - 1.5 * log(t) - LOG_2PI_2 - 0.5 * log(1 + sv*sv * t)
      + (sv*sv * a*a * w*w - 2 * v * a * w - v*v * t)
      / (2 + 2 * sv*sv * t);
    }
    return mult + log(sumf(t, a, w, 0, err / exp(mult)));
  }
}


// Gondan et al and Navarro et al
double nfb(const double& t, const double& a, const double& v,
          const double& w, const double& sv, const double& err,
          const double& sl_thresh, const NumFunc& numf, const SumFunc& sumf)
{ // note: sl_thresh is not used
  double mult;

  if (t / (a*a) > sl_thresh) { // use large-time
    if (sv <= SV_THRESH) { // no sv
      mult = exp(-v * a * w - v*v * t / 2) / (a*a);
    } else { // sv
      mult = exp((sv*sv * a*a * w*w - 2 * v * a * w - v*v * t)
                  / (2 + 2 * sv*sv * t)) / (a*a * sqrt(1 + sv*sv * t));
    }
    int kl = nkl_Nav(t / (a*a), w, err / mult);
    return mult * PI_CONST * nlarge_sum_Nav(t, a, w, kl, 0.0);
  } else { // use small-time
    if (sv <= SV_THRESH) { // no sv
      mult = exp(-v * a * w - v * v * t / 2);
    } else { // sv
      mult = exp((sv*sv * a*a * w*w - 2 * v * a * w - v*v * t)
                 / (2 + 2 * sv*sv * t)) / sqrt(1 + sv*sv * t);
    }
    int ks = numf(t / (a*a), w, err * a*a / mult);
    return mult * a * sumf(t, a, w, ks, 0.0)  / (t * SQRT_2PI * sqrt(t));
  }
}

double nfb_log(const double& t, const double& a, const double& v,
              const double& w, const double& sv, const double& err,
              const double& sl_thresh,
              const NumFunc& numf, const SumFunc& sumf)
{ // note: sl_thresh is not used
  double mult;

  if (t / (a*a) > sl_thresh) { // use large time
    if (sv <= SV_THRESH) { // no sv
      mult = - v * a * w - v * v * t / 2 - 2 * log(a);
    } else { // sv
      mult = (sv*sv * a*a * w*w - 2 * v * a * w - v*v * t)
      / (2 + 2 * sv*sv * t) - 0.5 * log(1 + sv*sv * t) - 2 * log(a);
    }
    int kl = nkl_Nav(t / (a*a), w, err / exp(mult));
    return LOG_PI + mult + log(nlarge_sum_Nav(t, a, w, kl, 0.0));
  } else { // use small time
    if (sv <= SV_THRESH) { // no sv
      mult = -v * a * w - v * v * t / 2;
    } else { // sv
      mult = (sv*sv * a*a * w*w - 2 * v * a * w - v*v * t)
             / (2 + 2 * sv*sv * t) - 0.5 * log(1 + sv*sv * t);
    }
    int ks = numf(t / (a*a), w, err * a*a / exp(mult));
    return mult + log(a) - LOG_2PI_2 - 1.5 * log(t)
           + log(sumf(t, a, w, ks, 0.0));
  }
}
