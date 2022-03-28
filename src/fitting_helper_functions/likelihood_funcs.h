// Functions to evaluate the DDM PDF and log(PDF)
// included from src/fitting_helper_functions/class_methods.h



// Density Function (non-log)
double ft(const double& t, const double& a, const double& v,
          const double& w, const double& sv, const double& err,
          const double& switch_thresh)
{
  double mult, sum_err;
  double taa = t / (a*a);

  if (taa > switch_thresh) { // use large time
    if (sv <= SV_THRESH) { // no sv
      mult = exp(-v * a * w - v*v * t / 2) / (a*a);
    } else { // sv
      mult = exp((sv*sv * a*a * w*w - 2 * v * a * w - v*v * t)
                  / (2 + 2 * sv*sv * t)) / (a*a * sqrt(1 + sv*sv * t));
    }
    // modify large-time error tolerance and check underflow
    sum_err = err / mult;
    if (sum_err < ERR_TOL_THRESH) sum_err = ERR_TOL_THRESH;
    int kl = kl_pdf(taa, sum_err);
    return PI_CONST * mult * large_sum(taa, w, kl);
  } else { // use small time
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
    return mult * small_sum(taa, w, sum_err);
  }
}

// Density Function (log)
double ft_log(const double& t, const double& a, const double& v,
              const double& w, const double& sv, const double& err,
              const double& switch_thresh)
{
  double mult, sum_err;
  double taa = t / (a*a);

  if (taa > switch_thresh) { // use large time
    if (sv <= SV_THRESH) { // no sv
      mult = - v * a * w - v * v * t / 2 - 2 * log(a);
    } else { // sv
      mult = (sv*sv * a*a * w*w - 2 * v * a * w - v*v * t)
      / (2 + 2 * sv*sv * t) - 0.5 * log(1 + sv*sv * t) - 2 * log(a);
    }
    // modify large-time error tolerance and check for underflow
    sum_err = err / exp(mult);
    if (sum_err < ERR_TOL_THRESH) sum_err = ERR_TOL_THRESH;
    int kl = kl_pdf(taa, sum_err);
    return LOG_PI + mult + log(large_sum(taa, w, kl));
  } else { // use small time
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
    return mult + log(small_sum(taa, w, sum_err));
  }
}
