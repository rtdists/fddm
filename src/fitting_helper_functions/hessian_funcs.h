// Functions to evaluate the Hessian of the DDM PDF
// included from src/fitting_helper_functions/class_methods.h



double dv2(const double& t, const double& v, const double& a, const double& w,
           const double& sv, const double& err, const double& sl_thresh)
{
  double taa = t / (a*a);
  double nnt = (1 + sv*sv*t);
  double onnt = 1 / nnt;
  double sqtonnt = sqrt(onnt);
  double mexp = exp(0.5 * onnt * (sv*sv*a*a*w*w - 2*v*a*w - v*v*t));
  double mult = mexp * onnt*onnt * sqtonnt * ((a*w + v*t)*(a*w + v*t) - t * nnt);
  double sum_err;

  if (taa > sl_thresh) { // use large-time
    mult /= a*a;
    sum_err = err / fabs(mult);
    if (sum_err < ERR_TOL_THRESH) sum_err = ERR_TOL_THRESH;
    int kl = kl_pdf(taa, sum_err);
    return mult * PI_CONST * large_sum(taa, w, kl);
  } else { // use small-time
    mult *= SQRT_1_2PI / (t * sqrt(taa));
    sum_err = err / fabs(mult);
    if (sum_err < ERR_TOL_THRESH) sum_err = ERR_TOL_THRESH;
    return mult * small_sum(taa, w, sum_err);
  }
}

double da2(const double& t, const double& v, const double& a, const double& w,
           const double& sv, const double& err, const double& sl_thresh)
{

  double taa = t / (a*a);
  double ot = 1 / t;
  double sqtot = sqrt(ot);
  double nnt = 1 + sv*sv * t;
  double onnt = 1 / nnt;
  double sqtonnt = sqrt(onnt);
  double mexp = exp(0.5 * onnt * (sv*sv*a*a*w*w - 2*v*a*w - v*v*t));
  double nawv = sv*sv*a*w - v;
  double arg = a*w*nawv + nnt;
  double gamma = sv*sv*a*a*w*w - v*a*w;
  double m0, m1, m2, sum_err0, sum_err1, sum_err2;

  if (taa > sl_thresh) { // use large-time
    m0 = mexp * onnt*onnt*sqtonnt / (a*a*a*a) * 
         (nnt*sv*sv*a*a*w*w + gamma*gamma - 4*gamma*nnt + 6*nnt*nnt);
    m1 = mexp * onnt*sqtonnt / (a*a*a) * (2*gamma - 7*nnt);
    m2 = mexp * sqtonnt * PI5 * taa*taa / (a*a*a*a);
    sum_err0 = err / fabs(m0);
    if (sum_err0 < ERR_TOL_THRESH) sum_err0 = ERR_TOL_THRESH;
    sum_err1 = err / fabs(m1);
    if (sum_err1 < ERR_TOL_THRESH) sum_err1 = ERR_TOL_THRESH;
    sum_err2 = err / fabs(m2);
    if (sum_err2 < ERR_TOL_THRESH) sum_err2 = ERR_TOL_THRESH;
    int kl0 = kl_pdf(taa, 0.33 * sum_err0);
    int kl1 = kl_dat(taa, t, 0.33 * sum_err1);
    int kl2 = kl_dat2(taa, 0.33 * sum_err2);
    return m0 * PI_CONST * large_sum(taa, w, kl0) +
           m1 * PI3 * taa/a * large_sum_dat(taa, w, kl1)
           + m2 * large_sum_dat2(taa, w, kl2);
  } else { // use small-time
    m0 = mexp * SQRT_1_2PI * ot*sqtot * onnt*onnt*sqtonnt * w *
         (nnt * (2*sv*sv*a*w - v) + nawv * arg);
    m1 = -mexp * SQRT_1_2PI * ot*ot*sqtot * onnt*sqtonnt * a * (2*arg + nnt);
    m2 = mexp * SQRT_1_2PI * ot*ot*ot*sqtot * sqtonnt * a*a*a;
    sum_err0 = err / fabs(m0);
    if (sum_err0 < ERR_TOL_THRESH) sum_err0 = ERR_TOL_THRESH;
    sum_err1 = err / fabs(m1);
    if (sum_err1 < ERR_TOL_THRESH) sum_err1 = ERR_TOL_THRESH;
    sum_err2 = err / fabs(m2);
    if (sum_err2 < ERR_TOL_THRESH) sum_err2 = ERR_TOL_THRESH;
    return m0 * small_sum(taa, w, 0.33 * sum_err0) +
           m1 * small_sum_dat(taa, w, 0.33 * sum_err1) +
           m2 * small_sum_dat2(taa, w, 0.33 * sum_err2);
  }
}

double dt02(const double& t, const double& v, const double& a, const double& w,
            const double& sv, const double& err, const double& sl_thresh)
{ // small-time is more stable
  double taa = t / (a*a);
  double ot = 1 / t;
  double sqtot = sqrt(ot);
  double nnt = 1 + sv*sv * t;
  double nnt34 = 3 + 4 * sv*sv * t;
  double nnt35 = 3 + 5 * sv*sv * t;
  double nnt38 = 3 + 8 * sv*sv * t;
  double nnt58 = 5 + 8 * sv*sv * t;
  double nnt78 = 7 + 8 * sv*sv * t;
  double onnt = 1 / nnt;
  double sqtonnt = sqrt(onnt);
  double nawv = sv*sv*a*w - v;
  double nawvaw = sv*sv*a*a*w*w - v*a*w;
  double naw2vaw = sv*sv*a*a*w*w - 2*v*a*w;
  double arg = naw2vaw - v*v*t;
  double phi = (sv*sv*a*w - v)*(sv*sv*a*w - v);
  double phit = phi * t;
  double psi = phi + nnt;
  double beta = nnt*nnt34 + phit;
  double mexp = exp(0.5 * onnt * arg);
  double m0, m1, m2, sum_err0, sum_err1, sum_err2;
  
  if (taa > sl_thresh) { // use large-time
    m0 = -0.25 * mexp * onnt*onnt*onnt*onnt*sqtonnt / (a*a) *
         (2*sv*sv*sv*sv*nnt*nnt - 5*sv*sv*nnt*psi - phi*psi);
    m1 = -mexp * onnt*onnt*sqtonnt * (sv*sv*nnt + phi) / (a*a);
    m2 = 0.25 * mexp * PI5 * sqtonnt / (a*a*a*a*a*a);
    sum_err0 = err / fabs(m0);
    if (sum_err0 < ERR_TOL_THRESH) sum_err0 = ERR_TOL_THRESH;
    sum_err1 = err / fabs(m1);
    if (sum_err1 < ERR_TOL_THRESH) sum_err1 = ERR_TOL_THRESH;
    sum_err2 = err / fabs(m2);
    if (sum_err2 < ERR_TOL_THRESH) sum_err2 = ERR_TOL_THRESH;
    int kl0 = kl_pdf(taa, 0.33 * sum_err0);
    int kl1 = kl_dat(taa, t, 0.33 * sum_err1);
    int kl2 = kl_dat2(taa, 0.33 * sum_err2);
    return m0 * PI_CONST * large_sum(taa, w, kl0) -
           m1 * 0.5 * PI3 / (a*a) * large_sum_dat(taa, w, kl1) +
           m2 * large_sum_dat2(taa, w, kl2);
  } else { // use small-time
    m0 = 0.25 * mexp * SQRT_1_2PI * a * ot*ot*ot*sqtot *
         onnt*onnt*onnt*onnt*sqtonnt *
         (phit*beta + 5*nnt*nnt*beta + 5*sv*sv*t*nnt*beta - 2*nnt*nnt*phit -
          2*sv*sv*t*nnt*nnt*nnt78);
    m1 = -0.25 * mexp * SQRT_1_2PI * a*a*a * ot*ot*ot*ot*sqtot *
         onnt*onnt*sqtonnt *
         (7*nnt*nnt + nnt*nnt35 + 2*v*v*t + 2*sv*sv*t*naw2vaw);
    m2 = 0.25 * mexp * SQRT_1_2PI * a*a*a*a*a * ot*ot*ot*ot*ot*sqtot * sqtonnt;
    sum_err0 = err / fabs(m0);
    if (sum_err0 < ERR_TOL_THRESH) sum_err0 = ERR_TOL_THRESH;
    sum_err1 = err / fabs(m1);
    if (sum_err1 < ERR_TOL_THRESH) sum_err1 = ERR_TOL_THRESH;
    sum_err2 = err / fabs(m2);
    if (sum_err2 < ERR_TOL_THRESH) sum_err2 = ERR_TOL_THRESH;
    return m0 * small_sum(taa, w, 0.33 * sum_err0) +
           m1 * small_sum_dat(taa, w, 0.33 * sum_err1) +
           m2 * small_sum_dat2(taa, w, 0.33 * sum_err2);
  }
}

double dw2(const double& t, const double& v, const double& a, const double& w,
           const double& sv, const double& err, const double& sl_thresh)
{
  double taa = t / (a*a);
  double ot = 1 / t;
  double sqtot = sqrt(ot);
  double nnt = 1 + sv*sv * t;
  double onnt = 1 / nnt;
  double sqtonnt = sqrt(onnt);
  double mexp = exp(0.5 * onnt * (sv*sv*a*a*w*w - 2*v*a*w - v*v*t));
  double nawv = sv*sv*a*w - v;
  double m0, m1, m2, sum_err0, sum_err1, sum_err2;
  
  if (taa > sl_thresh) { // use large-time
    m0 = mexp * onnt*onnt*sqtonnt * (sv*sv*nnt + nawv*nawv);
    m1 = 2 * mexp * onnt*sqtonnt * nawv / a;
    m2 = -mexp * sqtonnt * a / t;
    sum_err0 = err / fabs(m0);
    if (sum_err0 < ERR_TOL_THRESH) sum_err0 = ERR_TOL_THRESH;
    int kl0 = kl_pdf(taa, 0.33 * sum_err0);
    sum_err1 = err / fabs(m1);
    if (sum_err1 < ERR_TOL_THRESH) sum_err1 = ERR_TOL_THRESH;
    int kl1 = kl_dw(taa, w, 0.33 * sum_err1);
    sum_err2 = err / fabs(m2);
    if (sum_err2 < ERR_TOL_THRESH) sum_err2 = ERR_TOL_THRESH;
    int kl2 = kl_dat(taa, t, 0.33 * sum_err2);
    return m0 * PI_CONST * large_sum(taa, w, kl0) +
           m1 * PI2 * large_sum_dw(taa, w, kl1) +
           m2 * PI3 * taa / a * large_sum_dat(taa, w, kl2);
  } else { // use small-time
    m0 = mexp * SQRT_1_2PI * ot*ot*sqtot * onnt*onnt*sqtonnt * a*a*a *
         (t*nawv*nawv + sv*sv*t*nnt - 3*nnt*nnt);
    m1 = 2 * mexp * onnt*sqtonnt * nawv / a;
    m2 = mexp * SQRT_1_2PI * ot*ot*ot*sqtot * sqtonnt * a*a*a*a*a;
    sum_err0 = err / fabs(m0);
    if (sum_err0 < ERR_TOL_THRESH) sum_err0 = ERR_TOL_THRESH;
    sum_err1 = err / fabs(m1);
    if (sum_err1 < ERR_TOL_THRESH) sum_err1 = ERR_TOL_THRESH;
    int ks1 = ks_dw(taa, w, 0.33 * sum_err1);
    sum_err2 = err / fabs(m2);
    if (sum_err2 < ERR_TOL_THRESH) sum_err2 = ERR_TOL_THRESH;
    return m0 * small_sum(taa, w, 0.33 * sum_err0) +
           m1 * SQRT_1_2PI *ot*sqtot * a*a*a * small_sum_dw(taa, w, ks1) +
           m2 * small_sum_dat(taa, w, 0.33 * sum_err2);
  }
}

double dsv2(const double& t, const double& v, const double& a, const double& w,
            const double& sv, const double& err, const double& sl_thresh)
{
  double taa = t / (a*a);
  double nnt = (1 + sv*sv * t);
  double onnt = 1 / nnt;
  double sqtonnt = sqrt(onnt);
  double arg = sv*sv*a*a*w*w - 2*v*a*w - v*v*t;
  double mexp = exp(0.5 * onnt * arg);
  double gamma = a*a*w*w + 2*v*a*w*t + v*v*t*t - t*nnt;
  double mult = mexp * onnt*onnt*sqtonnt *
                (gamma - 2*sv*sv*t*t - 5*sv*sv*t*gamma * onnt +
                 sv*sv*gamma * onnt * (a*a*w*w - t*arg*onnt));
  double sum_err;

  if (taa > sl_thresh) { // use large-time
    mult /= a*a;
    sum_err = err / fabs(mult);
    if (sum_err < ERR_TOL_THRESH) sum_err = ERR_TOL_THRESH;
    int kl = kl_pdf(taa, sum_err);
    return mult * PI_CONST * large_sum(taa, w, kl);
  } else { // use small-time
    mult *= SQRT_1_2PI / (t * sqrt(taa));
    sum_err = err / fabs(mult);
    if (sum_err < ERR_TOL_THRESH) sum_err = ERR_TOL_THRESH;
    return mult * small_sum(taa, w, sum_err);
  }
}
