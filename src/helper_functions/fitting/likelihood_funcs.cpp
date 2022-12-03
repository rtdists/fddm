// Functions to evaluate the DDM PDF

#include "declarations.h"


// Density Function (non-log)
double pdf(const double& t, const double& v, const double& a, const double& w,
           const double& sv, const double& err, const double& switch_thresh)
{
  double mult, sum_err;
  double taa = t / (a*a);

  if (taa > switch_thresh) { // use large time
    mult = exp((sv*sv * a*a * w*w - 2 * v * a * w - v*v * t)
               / (2 + 2 * sv*sv * t)) / (a*a * sqrt(1 + sv*sv * t));
    // modify large-time error tolerance and check underflow
    sum_err = err / mult;
    if (sum_err < ERR_TOL_THRESH) sum_err = ERR_TOL_THRESH;
    int kl = kl_pdf(taa, sum_err);
    return PI_CONST * mult * large_sum(taa, w, kl);
  } else { // use small time
    mult = a * exp((sv*sv * a*a * w*w - 2 * v * a * w - v*v * t)
                   / (2 + 2 * sv*sv * t))
           / (t * SQRT_2PI * sqrt(t + sv*sv * t*t));
    // modify small-time error tolerance and check for underflow
    sum_err = err / mult;
    if (sum_err < ERR_TOL_THRESH) sum_err = ERR_TOL_THRESH;
    return mult * small_sum(taa, w, sum_err);
  }
}
