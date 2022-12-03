// Functions to evaluate the partial derivatives of the DDM PDF

#include "declarations.h"


double dv(const double& t, const double& v, const double& a, const double& w,
          const double& sv, const double& err, const double& sl_thresh)
{
  double taa = t / (a*a);
  double nnt = 1 / (1 + sv*sv * t);
  double sqtnnt = sqrt(nnt);
  double arg = a * w + v * t;
  double mexp = exp(0.5 * (sv*sv * a*a * w*w - 2 * v * a * w - v*v * t) * nnt);
  double mult, sum_err;

  if (taa > sl_thresh) { // use large-time
    mult = -mexp * arg * nnt * sqtnnt / (a*a);
    sum_err = err / fabs(mult);
    if (sum_err < ERR_TOL_THRESH) sum_err = ERR_TOL_THRESH;
    int kl = kl_pdf(taa, sum_err);
    return mult * PI_CONST * large_sum(taa, w, kl);
  } else { // use small-time
    mult = -mexp * arg * a * SQRT_1_2PI * nnt * sqtnnt / (t * sqrt(t));
    sum_err = err / fabs(mult);
    if (sum_err < ERR_TOL_THRESH) sum_err = ERR_TOL_THRESH;
    return mult * small_sum(taa, w, sum_err);
  }
}

double da(const double& t, const double& v, const double& a, const double& w,
          const double& sv, const double& err, const double& sl_thresh)
{
  double taa = t / (a*a);
  double nnt = 1 + sv*sv * t;
  double sqtnnt = sqrt(nnt);
  double naw = sv*sv * a*a * w*w;
  double vaw = v * a * w;
  double vvt = v*v * t;
  double mexp = exp(0.5 * (naw - 2 * vaw - vvt) / nnt);
  double m1, m2, sum_err1, sum_err2;

  if (taa > sl_thresh) { // use large-time
    m1 = mexp * (naw - vaw - 2 * nnt) / (a*a*a * nnt * sqtnnt);
    m2 = mexp / (a*a * sqtnnt);
    sum_err1 = err / fabs(m1);
    if (sum_err1 < ERR_TOL_THRESH) sum_err1 = ERR_TOL_THRESH;
    sum_err2 = err / fabs(m2);
    if (sum_err2 < ERR_TOL_THRESH) sum_err2 = ERR_TOL_THRESH;
    int kl1 = kl_pdf(t / (a*a), 0.5 * sum_err1);
    int kl2 = kl_dat(t / (a*a), t, 0.5 * sum_err2);
    return m1 * PI_CONST * large_sum(taa, w, kl1) +
           m2 * PI_CONST*PI_CONST*PI_CONST * t / (a*a*a) *
           large_sum_dat(taa, w, kl2);
  } else { // use small-time
    m1 = mexp * (naw - vaw + nnt) * SQRT_1_2PI / (t * sqrt(t) * nnt * sqtnnt);
    m2 = -mexp * a*a * SQRT_1_2PI / (t*t * sqrt(t) * sqtnnt);
    sum_err1 = err / fabs(m1);
    if (sum_err1 < ERR_TOL_THRESH) sum_err1 = ERR_TOL_THRESH;
    sum_err2 = err / fabs(m2);
    if (sum_err2 < ERR_TOL_THRESH) sum_err2 = ERR_TOL_THRESH;
    return m1 * small_sum(taa, w, 0.5 * sum_err1) +
           m2 * small_sum_dat(taa, w, 0.5 * sum_err2);
  }
}

double dt(const double& t, const double& v, const double& a, const double& w,
           const double& sv, const double& err, const double& sl_thresh)
{ // negative of dt0
  double taa = t / (a*a);
  double nnt = 1 + sv*sv * t;
  double sqtnnt = sqrt(nnt);
  double nawvaw = sv*sv * a*a * w*w - 2 * v * a * w;
  double mexp = exp(0.5 * (nawvaw - v*v * t) / nnt);
  double m1, m2, sum_err1, sum_err2;

  if (taa > sl_thresh) { // use large-time
    m1 = -0.5 * mexp * (v*v + sv*sv * (nawvaw + nnt)) / (a*a * nnt*nnt * sqtnnt);
    m2 = mexp / (a*a * sqtnnt);
    sum_err1 = err / fabs(m1);
    if (sum_err1 < ERR_TOL_THRESH) sum_err1 = ERR_TOL_THRESH;
    sum_err2 = err / fabs(m2);
    if (sum_err2 < ERR_TOL_THRESH) sum_err2 = ERR_TOL_THRESH;
    int kl1 = kl_pdf(t / (a*a), 0.5 * sum_err1);
    int kl2 = kl_dat(t / (a*a), t, 0.5 * sum_err2);
    return m1 * PI_CONST * large_sum(taa, w, kl1) -
           0.5 * m2 * PI_CONST*PI_CONST*PI_CONST / (a*a) *
           large_sum_dat(taa, w, kl2);
  } else { // use small-time
    m1 = -0.5 * mexp * SQRT_1_2PI * a *
         (nnt * (3 + 4 * sv*sv * t) + v*v * t + sv*sv*t * nawvaw) /
         (t*t * sqrt(t) * nnt*nnt * sqtnnt);
    m2 = 0.5 * mexp * SQRT_1_2PI * a*a*a / (t*t*t * sqrt(t) * sqtnnt);
    sum_err1 = err / fabs(m1);
    if (sum_err1 < ERR_TOL_THRESH) sum_err1 = ERR_TOL_THRESH;
    sum_err2 = err / fabs(m2);
    if (sum_err2 < ERR_TOL_THRESH) sum_err2 = ERR_TOL_THRESH;
    return m1 * small_sum(taa, w, 0.5 * sum_err1) +
           m2 * small_sum_dat(taa, w, 0.5 * sum_err2);
  }
}

double dt0(const double& t, const double& v, const double& a, const double& w,
           const double& sv, const double& err, const double& sl_thresh)
{ // negative of dt
  double taa = t / (a*a);
  double nnt = 1 + sv*sv * t;
  double sqtnnt = sqrt(nnt);
  double nawvaw = sv*sv * a*a * w*w - 2 * v * a * w;
  double mexp = exp(0.5 * (nawvaw - v*v * t) / nnt);
  double m1, m2, sum_err1, sum_err2;

  if (taa > sl_thresh) { // use large-time
    m1 = 0.5 * mexp * (v*v + sv*sv * (nawvaw + nnt)) / (a*a * nnt*nnt * sqtnnt);
    m2 = -mexp / (a*a * sqtnnt);
    sum_err1 = err / fabs(m1);
    if (sum_err1 < ERR_TOL_THRESH) sum_err1 = ERR_TOL_THRESH;
    sum_err2 = err / fabs(m2);
    if (sum_err2 < ERR_TOL_THRESH) sum_err2 = ERR_TOL_THRESH;
    int kl1 = kl_pdf(t / (a*a), 0.5 * sum_err1);
    int kl2 = kl_dat(t / (a*a), t, 0.5 * sum_err2);
    return m1 * PI_CONST * large_sum(taa, w, kl1) -
           0.5 * m2 * PI_CONST*PI_CONST*PI_CONST / (a*a) *
           large_sum_dat(taa, w, kl2);
  } else { // use small-time
    m1 = 0.5 * mexp * SQRT_1_2PI * a *
         (nnt * (3 + 4 * sv*sv * t) + v*v * t + sv*sv*t * nawvaw) /
         (t*t * sqrt(t) * nnt*nnt * sqtnnt);
    m2 = -0.5 * mexp * SQRT_1_2PI * a*a*a / (t*t*t * sqrt(t) * sqtnnt);
    sum_err1 = err / fabs(m1);
    if (sum_err1 < ERR_TOL_THRESH) sum_err1 = ERR_TOL_THRESH;
    sum_err2 = err / fabs(m2);
    if (sum_err2 < ERR_TOL_THRESH) sum_err2 = ERR_TOL_THRESH;
    return m1 * small_sum(taa, w, 0.5 * sum_err1) +
           m2 * small_sum_dat(taa, w, 0.5 * sum_err2);
  }
}

double dw(const double& t, const double& v, const double& a, const double& w,
          const double& sv, const double& err, const double& sl_thresh)
{
  double out;
  double taa = t / (a*a);
  double onnt = 1 / (1 + sv*sv * t);
  double osqtnnt = sqrt(onnt);
  double arg = sv*sv * a * w - v;
  double mexp = exp(0.5 * onnt * (sv*sv * a*a * w*w - 2 * v * a * w - v*v * t));
  double mult, sum_err;

  // first part
  mult = mexp * onnt * osqtnnt * arg / a;
  sum_err = err / fabs(mult);
  if (sum_err < ERR_TOL_THRESH) sum_err = ERR_TOL_THRESH;
  int kl1 = kl_pdf(t / (a*a), 0.5 * sum_err);
  if (kl1 <= sl_thresh) { // use large-time (same heuristic from dfddm)
    out = mult * PI_CONST * large_sum(taa, w, kl1);
  } else { // use small-time
    mult = mexp * arg * a*a * SQRT_1_2PI * onnt * osqtnnt / (t * sqrt(t));
    out = mult * small_sum(taa, w, 0.5 * sum_err);
  }

  // second part
  mult = mexp * osqtnnt / (a*a);
  sum_err = err / fabs(mult);
  if (sum_err < ERR_TOL_THRESH) sum_err = ERR_TOL_THRESH;
  int kl2 = kl_dw(t / (a*a), t, 0.5 * sum_err);
  int ks2 = ks_dw(t / (a*a), w, 0.5 * sum_err);
  if (kl2 < 2 * ks2) { // use large-time
    out += mult * PI_CONST*PI_CONST * large_sum_dw(taa, w, kl2);
  } else { // use small time
    out += mult * a*a*a * SQRT_1_2PI / (t * sqrt(t)) * small_sum_dw(taa, w, ks2);
  }

  return out;
}

double dsv(const double& t, const double& v, const double& a, const double& w,
           const double& sv, const double& err, const double& sl_thresh)
{
  double taa = t / (a*a);
  double nnt = 1 + sv*sv * t;
  double sqtnnt = sqrt(nnt);
  double arg = a*a * w*w + 2 * v * a * w * t + v*v * t*t - t * nnt;
  double mexp = exp(0.5 * (sv*sv * a*a * w*w - 2 * v * a * w - v*v * t) / nnt);
  double mult, sum_err;

  if (taa > sl_thresh) { // use large-time
    mult = mexp * arg * sv / (a*a * nnt*nnt * sqtnnt);
    sum_err = err / fabs(mult);
    if (sum_err < ERR_TOL_THRESH) sum_err = ERR_TOL_THRESH;
    int kl = kl_pdf(taa, sum_err);
    return mult * PI_CONST * large_sum(taa, w, kl);
  } else { // use small-time
    mult = mexp * arg * sv * a * SQRT_1_2PI / (t * sqrt(t) * nnt*nnt * sqtnnt);
    sum_err = err / fabs(mult);
    if (sum_err < ERR_TOL_THRESH) sum_err = ERR_TOL_THRESH;
    return mult * small_sum(taa, w, sum_err);
  }
}
