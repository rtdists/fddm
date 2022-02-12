// Functions to evaluate the DDM PDF for specific criteria



//////////                                                            //////////
///////////////////////////////// Small Time ///////////////////////////////////
//////////                                                            //////////


// Stop When Small Enough
double ff(const double& t, const double& a, const double& v,
          const double& w, const double& sv, const double& err,
          const int& max_terms_large, const NumFunc& numf, const SumFunc& sumf)
{ // note: max_terms_large, numf are not used
  double mult, sum_err;

  // calculate multiplicative term outside sum
  if (sv <= SV_THRESH) { // no sv
    mult = a * exp(-v * a * w - v * v * t / 2) / (t * SQRT_2PI * sqrt(t));
  } else { // sv
    mult = a * exp((sv*sv * a*a * w*w - 2 * v * a * w - v*v * t)
                     / (2 + 2 * sv*sv * t)) / (t * SQRT_2PI * sqrt(t + sv*sv * t*t));
  }

  // modify error tolerance and check for underflow
  sum_err = err / mult;
  if (sum_err < ERR_TOL_THRESH) sum_err = ERR_TOL_THRESH;
  return mult * sumf(t, a, w, 0, sum_err);
}

double ff_log(const double& t, const double& a, const double& v,
              const double& w, const double& sv, const double& err,
              const int& max_terms_large,
              const NumFunc& numf, const SumFunc& sumf)
{ // note: max_terms_large, numf are not used
  double mult, sum_err;

  // calculate multiplicative term outside sum
  if (sv <= SV_THRESH) { // no sv
    mult = log(a) - LOG_2PI_2 - 1.5 * log(t) - v * a * w - v*v * t / 2;
  } else { // sv
    mult = log(a) - 1.5 * log(t) - LOG_2PI_2 - 0.5 * log(1 + sv*sv * t)
    + (sv*sv * a*a * w*w - 2 * v * a * w - v*v * t)
    / (2 + 2 * sv*sv * t);
  }

  // modify error tolerance and check for underflow
  sum_err = err / exp(mult);
  if (sum_err < ERR_TOL_THRESH) sum_err = ERR_TOL_THRESH;
  return mult + log(sumf(t, a, w, 0, sum_err));
}


// Gondan et al and Navarro et al
double fs(const double& t, const double& a, const double& v,
          const double& w, const double& sv, const double& err,
          const int& max_terms_large, const NumFunc& numf, const SumFunc& sumf)
{ // note: max_terms_large is not used
  double mult_s, sum_err;
  int ks;

  // calculate multiplicative term outside sum
  if (sv <= SV_THRESH) { // no sv
    mult_s = exp(-v * a * w - v * v * t / 2);
  } else { // sv
    mult_s = exp((sv*sv * a*a * w*w - 2 * v * a * w - v*v * t)
                   / (2 + 2 * sv*sv * t)) / sqrt(1 + sv*sv * t);
  }

  // modify error tolerance and check for underflow
  sum_err = err / mult_s;
  if (sum_err < ERR_TOL_THRESH) sum_err = ERR_TOL_THRESH;
  ks = numf(t / (a*a), w, sum_err * a*a);
  return mult_s * a * sumf(t, a, w, ks, 0.0) / (t * SQRT_2PI * sqrt(t));
}

double fs_log(const double& t, const double& a, const double& v,
              const double& w, const double& sv, const double& err,
              const int& max_terms_large,
              const NumFunc& numf, const SumFunc& sumf)
{ // note: max_terms_large is not used
  double mult_s, sum_err;
  int ks;

  // calculate multiplicative term outside sum
  if (sv <= SV_THRESH) { // no sv
    mult_s = -v * a * w - v * v * t / 2;
  } else { // sv
    mult_s = (sv*sv * a*a * w*w - 2 * v * a * w - v*v * t)
    / (2 + 2 * sv*sv * t) - 0.5 * log(1 + sv*sv * t);
  }

  // modify error tolerance and check for underflow
  sum_err = err / exp(mult_s);
  if (sum_err < ERR_TOL_THRESH) sum_err = ERR_TOL_THRESH;
  ks = numf(t / (a*a), w, sum_err * a*a);
  return mult_s + log(a) + log(sumf(t, a, w, ks, 0.0)) - LOG_2PI_2
         - 1.5 * log(t);
}





//////////                                                            //////////
///////////////////////////////// Large Time ///////////////////////////////////
//////////                                                            //////////


// Navarro et al
double fl(const double& t, const double& a, const double& v,
          const double& w, const double& sv, const double& err,
          const int& max_terms_large, const NumFunc& numf, const SumFunc& sumf)
{ // note: max_terms_large, numf, sumf are not used
  double mult_l, sum_err;
  int kl;

  // calculate multiplicative term outside sum
  if (sv <= SV_THRESH) { // no sv
    mult_l = exp(-v * a * w - v*v * t / 2) / (a*a);
  } else { // sv
    mult_l = exp((sv*sv * a*a * w*w - 2 * v * a * w - v*v * t)
                   / (2 + 2 * sv*sv * t)) / (a*a * sqrt(1 + sv*sv * t));
  }

  // modify error tolerance and check for underflow
  sum_err = err / mult_l;
  if (sum_err < ERR_TOL_THRESH) sum_err = ERR_TOL_THRESH;
  kl = kl_Nav(t / (a*a), w, sum_err);
  return mult_l * PI_CONST * large_sum_Nav(t, a, w, kl, 0.0);
}

double fl_log(const double& t, const double& a, const double& v,
              const double& w, const double& sv, const double& err,
              const int& max_terms_large,
              const NumFunc& numf, const SumFunc& sumf)
{ // note: max_terms_large, numf, sumf are not used
  double mult_l, sum_err;
  int kl;

  // calculate multiplicative term outside sum
  if (sv <= SV_THRESH) { // no sv
    mult_l = - v * a * w - v * v * t / 2 - 2 * log(a);
  } else { // sv
    mult_l = (sv*sv * a*a * w*w - 2 * v * a * w - v*v * t)
    / (2 + 2 * sv*sv * t) - 0.5 * log(1 + sv*sv * t) - 2 * log(a);
  }

  // modify error tolerance and check for underflow
  sum_err = err / exp(mult_l);
  if (sum_err < ERR_TOL_THRESH) sum_err = ERR_TOL_THRESH;
  kl = kl_Nav(t / (a*a), w, sum_err);
  return mult_l + LOG_PI + log(large_sum_Nav(t, a, w, kl, 0.0));
}



//////////                                                            //////////
///////////////////////////////// Both Times ///////////////////////////////////
//////////                                                            //////////


// Stop When Small Enough
double fc(const double& t, const double& a, const double& v,
          const double& w, const double& sv, const double& err,
          const int& max_terms_large, const NumFunc& numf, const SumFunc& sumf)
{ // note: numf is not used
  double mult, sum_err;

  // calculate large-time multiplicative term outside sum
  if (sv <= SV_THRESH) { // no sv
    mult = exp(-v * a * w - v*v * t / 2) / (a*a);
  } else { // sv
    mult = exp((sv*sv * a*a * w*w - 2 * v * a * w - v*v * t)
                 / (2 + 2 * sv*sv * t)) / (a*a * sqrt(1 + sv*sv * t));
  }
  // modify large-time error tolerance and check underflow
  sum_err = err / mult;
  if (sum_err < ERR_TOL_THRESH) sum_err = ERR_TOL_THRESH;
  int kl = kl_Nav(t / (a*a), w, sum_err);

  // choose large vs small time
  if (kl <= max_terms_large) {
    return PI_CONST * mult * large_sum_Nav(t, a, w, kl, 0.0);
  } else {
    if (sv <= SV_THRESH) { // no sv
      mult = a * exp(-v * a * w - v * v * t / 2) / (t * SQRT_2PI * sqrt(t));
    } else { // sv
      mult = a * exp((sv*sv * a*a * w*w - 2 * v * a * w - v*v * t)
                       / (2 + 2 * sv*sv * t))
             / (t * SQRT_2PI * sqrt(t + sv*sv * t*t));
    }
    // modify small-time error tolerance and check for underflow
    sum_err = err / mult;
    if (sum_err < ERR_TOL_THRESH) sum_err = ERR_TOL_THRESH;
    return mult * sumf(t, a, w, 0, sum_err);
  }
}

double fc_log(const double& t, const double& a, const double& v,
              const double& w, const double& sv, const double& err,
              const int& max_terms_large,
              const NumFunc& numf, const SumFunc& sumf)
{ // note: numf is not used
  double mult, sum_err;

  // calculate large-time multiplicative term outside sum
  if (sv <= SV_THRESH) { // no sv
    mult = - v * a * w - v * v * t / 2 - 2 * log(a);
  } else { // sv
    mult = (sv*sv * a*a * w*w - 2 * v * a * w - v*v * t)
    / (2 + 2 * sv*sv * t) - 0.5 * log(1 + sv*sv * t) - 2 * log(a);
  }
  // modify large-time error tolerance and check for underflow
  sum_err = err / exp(mult);
  if (sum_err < ERR_TOL_THRESH) sum_err = ERR_TOL_THRESH;
  int kl = kl_Nav(t / (a*a), w, sum_err);

  // choose large vs small time
  if (kl <= max_terms_large) {
    return LOG_PI + mult + log(large_sum_Nav(t, a, w, kl, 0.0));
  } else {
    if (sv <= SV_THRESH) { // no sv
      mult = log(a) - LOG_2PI_2 - 1.5 * log(t) - v * a * w - v*v * t / 2;
    } else { // sv
      mult = log(a) - 1.5 * log(t) - LOG_2PI_2 - 0.5 * log(1 + sv*sv * t)
      + (sv*sv * a*a * w*w - 2 * v * a * w - v*v * t)
      / (2 + 2 * sv*sv * t);
    }
    // modify small-time error tolerance and check for underflow
    sum_err = err / exp(mult);
    if (sum_err < ERR_TOL_THRESH) sum_err = ERR_TOL_THRESH;
    return mult + log(sumf(t, a, w, 0, sum_err));
  }
}


// Gondan et al and Navarro et al
double fb(const double& t, const double& a, const double& v,
          const double& w, const double& sv, const double& err,
          const int& max_terms_large, const NumFunc& numf, const SumFunc& sumf)
{ // note: max_terms_large is not used
  double mult_s, mult_l, sum_err_s, sum_err_l;
  int ks, kl;

  // calculate multiplicative term outside sum
  if (sv <= SV_THRESH) { // no sv
    mult_s = exp(-v * a * w - v * v * t / 2);
    mult_l = exp(-v * a * w - v*v * t / 2) / (a*a);
  } else { // sv
    mult_s = exp((sv*sv * a*a * w*w - 2 * v * a * w - v*v * t)
                   / (2 + 2 * sv*sv * t)) / sqrt(1 + sv*sv * t);
    mult_l = exp((sv*sv * a*a * w*w - 2 * v * a * w - v*v * t)
                   / (2 + 2 * sv*sv * t)) / (a*a * sqrt(1 + sv*sv * t));
  }

  // modify error tolerance and check for underflow
  sum_err_s = err / mult_s;
  if (sum_err_s < ERR_TOL_THRESH) sum_err_s = ERR_TOL_THRESH;
  ks = numf(t / (a*a), w, sum_err_s * a*a);
  sum_err_l = err / mult_l;
  if (sum_err_l < ERR_TOL_THRESH) sum_err_l = ERR_TOL_THRESH;
  kl = kl_Nav(t / (a*a), w, sum_err_l);

  // choose small vs large time
  if (ks < kl) {
    return mult_s * a * sumf(t, a, w, ks, 0.0) / (t * SQRT_2PI * sqrt(t));
  } else {
    return mult_l * PI_CONST * large_sum_Nav(t, a, w, kl, 0.0);
  }
}

double fb_log(const double& t, const double& a, const double& v,
              const double& w, const double& sv, const double& err,
              const int& max_terms_large,
              const NumFunc& numf, const SumFunc& sumf)
{ // note: max_terms_large is not used
  double mult_s, mult_l, sum_err_s, sum_err_l;
  int ks, kl;

  // calculate multiplicative term outside sum
  if (sv <= SV_THRESH) { // no sv
    mult_s = -v * a * w - v * v * t / 2;
    mult_l = - v * a * w - v * v * t / 2 - 2 * log(a);
  } else { // sv
    mult_s = (sv*sv * a*a * w*w - 2 * v * a * w - v*v * t)
    / (2 + 2 * sv*sv * t) - 0.5 * log(1 + sv*sv * t);
    mult_l = (sv*sv * a*a * w*w - 2 * v * a * w - v*v * t)
      / (2 + 2 * sv*sv * t) - 0.5 * log(1 + sv*sv * t) - 2 * log(a);
  }

  // modify error tolerance and check for underflow
  sum_err_s = err / exp(mult_s);
  if (sum_err_s < ERR_TOL_THRESH) sum_err_s = ERR_TOL_THRESH;
  ks = numf(t / (a*a), w, sum_err_s * a*a);
  sum_err_l = err / exp(mult_l);
  if (sum_err_l < ERR_TOL_THRESH) sum_err_l = ERR_TOL_THRESH;
  kl = kl_Nav(t / (a*a), w, sum_err_l);

  // choose small vs large time
  if (ks < kl) {
    return mult_s + log(a) + log(sumf(t, a, w, ks, 0.0)) - LOG_2PI_2
           - 1.5 * log(t);
  } else {
    return mult_l + LOG_PI + log(large_sum_Nav(t, a, w, kl, 0.0));
  }
}
