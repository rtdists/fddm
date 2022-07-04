// Functions to evaluate the gradient of the DDM PDF
// included from src/fitting_helper_functions/class_methods.h



double dv2(const double& t, const double& a, const double& v, const double& w,
           const double& sv, const double& err, const double& sl_thresh)
{
  double eps = 0.5 * err;
  return (dv(t, a, v+eps, w, sv, err, sl_thresh) -
          dv(t, a, v-eps, w, sv, err, sl_thresh)) / (err); // err = 2 * eps

  // double taa = t / (a*a);
  // double nnt = (1 + sv*sv * t);
  // double onnt = 1 / nnt;
  // double sqtonnt = sqrt(onnt);
  // double mexp = exp(0.5 * onnt * (sv*sv*a*a*w*w - 2*v*a*w - v*v*t));
  // double mult = onnt*onnt * sqtonnt * ((a*w + v*t)*(a*w + v*t) - t* nnt);
  // double sum_err;

  // if (taa > sl_thresh) { // use large-time
  //   mult /= a*a;
  //   sum_err = err / fabs(mult);
  //   if (sum_err < ERR_TOL_THRESH) sum_err = ERR_TOL_THRESH;
  //   int kl = kl_pdf(taa, sum_err);
  //   return mult * PI_CONST * large_sum(taa, w, kl);
  // } else { // use small-time
  //   mult *= SQRT_1_2PI * a / (t * sqrt(t));
  //   sum_err = err / fabs(mult);
  //   if (sum_err < ERR_TOL_THRESH) sum_err = ERR_TOL_THRESH;
  //   return mult * small_sum(taa, w, sum_err);
  // }
}

double da2(const double& t, const double& a, const double& v, const double& w,
           const double& sv, const double& err, const double& sl_thresh)
{
  double eps = 0.5 * err;
  while (a - eps <= 0) {
    eps *= 0.5;
  }
  return (da(t, a+eps, v, w, sv, err, sl_thresh) -
          da(t, a-eps, v, w, sv, err, sl_thresh)) / (err); // err = 2 * eps

  // double taa = t / (a*a);
  // double ot = 1 / t;
  // double sqtot = sqrt(ot);
  // double nnt = 1 + sv*sv * t;
  // double onnt = 1 / nnt;
  // double sqtonnt = sqrt(onnt);
  // double mexp = exp(0.5 * onnt * (sv*sv*a*a*w*w - 2*v*a*w - v*v*t));
  // double nawv = sv*sv*a*w - v;
  // double arg = a*w*nawv + nnt;
  // double argl = a*w*nawv - 2*nnt;
  // double m0, m1, m2, sum_err0, sum_err1, sum_err2;

  // if (taa > sl_thresh) { // use large-time
  //   m0 = mexp * onnt*onnt*sqtonnt / (a*a*a*a) * 
  //        (a*w*(2*sv*sv*a*w - v)*nnt + argl*(a*w*nawv - 3*nnt));
  //   m1 = 2 * mexp * onnt*sqtonnt * argl / (a*a*a);
  //   m2 = mexp * PI_CONST*PI_CONST*PI_CONST*PI_CONST*PI_CONST * sqtonnt *
  //       taa*taa / (a*a*a*a);
  //   sum_err0 = err / fabs(m0);
  //   if (sum_err0 < ERR_TOL_THRESH) sum_err0 = ERR_TOL_THRESH;
  //   sum_err1 = err / fabs(m1);
  //   if (sum_err1 < ERR_TOL_THRESH) sum_err1 = ERR_TOL_THRESH;
  //   sum_err2 = err / fabs(m2);
  //   if (sum_err2 < ERR_TOL_THRESH) sum_err2 = ERR_TOL_THRESH;
  //   int kl0 = kl_pdf(taa, 0.33 * sum_err0);
  //   int kl1 = kl_dat(taa, t, 0.33 * sum_err1);
  //   int kl2 = kl_dat2(taa, 0.33 * sum_err2);
  //   return m0 * PI_CONST * large_sum(taa, w, kl0) +
  //          m1 * PI_CONST*PI_CONST*PI_CONST * taa/a * large_sum_dat(taa, w, kl1)
  //          + m2 * large_sum_dat2(taa, w, kl2);
  // } else { // use small-time
  //   m0 = mexp * SQRT_1_2PI * ot*sqtot * onnt*onnt*sqtonnt * w *
  //        (nnt * (2*sv*sv*a*w - v) + nawv * arg);
  //   m1 = -mexp * SQRT_1_2PI * ot*ot*sqtot * onnt*sqtonnt * a * (2*arg + nnt);
  //   m2 = mexp * SQRT_1_2PI * ot*ot*ot*sqtot * sqtonnt * a*a*a;
  //   sum_err0 = err / fabs(m0);
  //   if (sum_err0 < ERR_TOL_THRESH) sum_err0 = ERR_TOL_THRESH;
  //   sum_err1 = err / fabs(m1);
  //   if (sum_err1 < ERR_TOL_THRESH) sum_err1 = ERR_TOL_THRESH;
  //   sum_err2 = err / fabs(m2);
  //   if (sum_err2 < ERR_TOL_THRESH) sum_err2 = ERR_TOL_THRESH;
  //   return m0 * small_sum(taa, w, 0.33 * sum_err0) +
  //          m1 * small_sum_dat(taa, w, 0.33 * sum_err1) +
  //          m2 * small_sum_dat2(taa, w, 0.33 * sum_err2);
  // }
}

double dt02(const double& t, const double& a, const double& v, const double& w,
            const double& sv, const double& err, const double& sl_thresh)
{
  double eps = 0.5 * err;
  while (t - eps <= 0) {
    eps *= 0.5;
  }
  return (dt0(t+eps, a, v, w, sv, err, sl_thresh) -
          dt0(t-eps, a, v, w, sv, err, sl_thresh)) / (err); // err = 2 * eps

  // double taa = t / (a*a);
  // double ot = 1 / t;
  // double sqtot = sqrt(ot);
  // double nnt = 1 + sv*sv * t;
  // double nnt34 = 3 + 4 * sv*sv * t;
  // double onnt = 1 / nnt;
  // double sqtonnt = sqrt(onnt);
  // double nawvaw = sv*sv*a*a*w*w - v*a*w;
  // double naw2vaw = sv*sv*a*a*w*w - 2*v*a*w;
  // double arg = naw2vaw - v*v*t;
  // double mexp = exp(0.5 * onnt * arg);
  // double m0, m1, m2, sum_err0, sum_err1, sum_err2;
  
  // if (taa > sl_thresh) { // use large-time
  //   m0 = 0.25 * mexp * onnt*onnt*onnt*onnt*sqtonnt / (a*a) * 
  //        (v*v*v*v + sv*sv*v*v*nnt + 3*sv*sv*sv*sv*nnt +
  //         sv*sv*naw2vaw*(2*v*v + 6*sv*sv*nnt + sv*sv*naw2vaw));
  //   m1 = -mexp * onnt*onnt*sqtonnt * (sv*sv*naw2vaw + (sv*sv + v*v)*nnt) /
  //        (a*a);
  //   m2 = 0.25 * mexp * PI_CONST*PI_CONST*PI_CONST*PI_CONST*PI_CONST * sqtonnt /
  //        (a*a*a*a*a*a);
  //   sum_err0 = err / fabs(m0);
  //   if (sum_err0 < ERR_TOL_THRESH) sum_err0 = ERR_TOL_THRESH;
  //   sum_err1 = err / fabs(m1);
  //   if (sum_err1 < ERR_TOL_THRESH) sum_err1 = ERR_TOL_THRESH;
  //   sum_err2 = err / fabs(m2);
  //   if (sum_err2 < ERR_TOL_THRESH) sum_err2 = ERR_TOL_THRESH;
  //   int kl0 = kl_pdf(taa, 0.33 * sum_err0);
  //   int kl1 = kl_dat(taa, t, 0.33 * sum_err1);
  //   int kl2 = kl_dat2(taa, 0.33 * sum_err2);
  //   return m0 * PI_CONST * large_sum(taa, w, kl0) -
  //          m1 * 0.5 * PI_CONST*PI_CONST*PI_CONST / (a*a) *
  //               large_sum_dat(taa, w, kl1) +
  //          m2 * large_sum_dat2(taa, w, kl2);
  // } else { // use small-time
  //   m0 = -0.25 * mexp * SQRT_1_2PI * a * ot*ot*ot*sqtot *
  //        onnt*onnt*onnt*onnt*sqtonnt *
  //        (t*nnt*nnt *
  //           (14*sv*sv - 3*v*v - 3*sv*sv*nawvaw - nnt34*(5*sv*sv*sv*sv + v*v)) -
  //         t*t*nnt *
  //           (5*sv*sv*v*v + v*v*v*v + sv*sv*(5*sv*sv*nawvaw + v*v*naw2vaw)) -
  //         sv*sv*t*arg * (nnt*nnt34 + v*v*t + sv*sv*t*naw2vaw) +
  //         4*sv*sv*sv*sv*t*t*nnt*nnt - 5*nnt*nnt*nnt*nnt34);
  //   m1 = -0.25 * mexp * SQRT_1_2PI * a*a*a * ot*ot*ot*sqtot *
  //        onnt*onnt*sqtonnt *
  //        (nnt*nnt34 + v*v*t + sv*sv*t*naw2vaw + ot*ot*ot*ot*sqtot * sqtonnt *
  //         (nnt*(14 + sv*sv*t + v*t) + sv*sv*t * arg));
  //   m2 = 0.125 * mexp * SQRT_1_2PI * a*a*a*a*a*a*a *
  //        ot*ot*ot*ot*ot*ot*ot*sqtot * sqtonnt;
  //   sum_err0 = err / fabs(m0);
  //   if (sum_err0 < ERR_TOL_THRESH) sum_err0 = ERR_TOL_THRESH;
  //   sum_err1 = err / fabs(m1);
  //   if (sum_err1 < ERR_TOL_THRESH) sum_err1 = ERR_TOL_THRESH;
  //   sum_err2 = err / fabs(m2);
  //   if (sum_err2 < ERR_TOL_THRESH) sum_err2 = ERR_TOL_THRESH;
  //   return m0 * small_sum(taa, w, 0.33 * sum_err0) +
  //          m1 * small_sum_dat(taa, w, 0.33 * sum_err1) +
  //          m2 * small_sum_dat2(taa, w, 0.33 * sum_err2);
  // }
}

double dw2(const double& t, const double& a, const double& v, const double& w,
           const double& sv, const double& err, const double& sl_thresh)
{
  double eps = 0.5 * err;
  while (w - eps <= 0 || w + eps >= 1) {
    eps *= 0.5;
  }
  return (dw(t, a, v, w+eps, sv, err, sl_thresh) -
          dw(t, a, v, w-eps, sv, err, sl_thresh)) / (err); // err = 2 * eps

  // double taa = t / (a*a);
  // double ot = 1 / t;
  // double sqtot = sqrt(ot);
  // double nnt = 1 + sv*sv * t;
  // double onnt = 1 / nnt;
  // double sqtonnt = sqrt(onnt);
  // double mexp = exp(0.5 * onnt * (sv*sv*a*a*w*w - 2*v*a*w - v*v*t));
  // double nawv = sv*sv*a*w - v;
  // double m0, m1, m2, sum_err0, sum_err1, sum_err2;
  
  // if (taa > sl_thresh) { // use large-time
  //   m0 = mexp * onnt*onnt*sqtonnt * (sv*sv*nnt + nawv*nawv);
  //   m1 = 2 * mexp * onnt*sqtonnt * nawv / a;
  //   m2 = mexp * sqtonnt / (a*a);
  //   sum_err0 = err / fabs(m0);
  //   if (sum_err0 < ERR_TOL_THRESH) sum_err0 = ERR_TOL_THRESH;
  //   int kl0 = kl_pdf(taa, 0.33 * sum_err0);
  //   sum_err1 = err / fabs(m1);
  //   if (sum_err1 < ERR_TOL_THRESH) sum_err1 = ERR_TOL_THRESH;
  //   int kl1 = kl_dw(taa, w, 0.33 * sum_err1);
  //   sum_err2 = err / fabs(m2);
  //   if (sum_err2 < ERR_TOL_THRESH) sum_err2 = ERR_TOL_THRESH;
  //   int kl2 = kl_dat(taa, t, 0.33 * sum_err2);
  //   return m0 * PI_CONST * large_sum(taa, w, kl0) +
  //          m1 * PI_CONST*PI_CONST * large_sum_dw(taa, w, kl1) -
  //          m2 * a / taa * large_sum_dat(taa, w, kl2);
  // } else { // use small-time
  //   m0 = mexp * SQRT_1_2PI * ot*ot*sqtot * onnt*onnt*sqtonnt * a*a*a *
  //        (t*nawv*nawv + sv*sv*t*nnt - 3*nnt*nnt);
  //   m1 = 2 * mexp * onnt*sqtonnt * nawv / a;
  //   m2 = mexp * SQRT_1_2PI * ot*ot*ot*sqtot * sqtonnt * a*a*a*a*a;
  //   sum_err0 = err / fabs(m0);
  //   if (sum_err0 < ERR_TOL_THRESH) sum_err0 = ERR_TOL_THRESH;
  //   sum_err1 = err / fabs(m1);
  //   if (sum_err1 < ERR_TOL_THRESH) sum_err1 = ERR_TOL_THRESH;
  //   int ks1 = ks_dw(taa, w, 0.33 * sum_err1);
  //   sum_err2 = err / fabs(m2);
  //   if (sum_err2 < ERR_TOL_THRESH) sum_err2 = ERR_TOL_THRESH;
  //   return m0 * small_sum(taa, w, 0.33 * sum_err0) +
  //          m1 * SQRT_1_2PI *ot*sqtot * a*a*a * small_sum_dw(taa, w, ks1) +
  //          m2 * small_sum_dat(taa, w, 0.33 * sum_err2);
  // }
}

double dsv2(const double& t, const double& a, const double& v, const double& w,
            const double& sv, const double& err, const double& sl_thresh)
{
  double eps = 0.5 * err;
  while (sv - eps <= 0) {
    eps *= 0.5;
  }
  return (dsv(t, a, v, w, sv+eps, err, sl_thresh) -
          dsv(t, a, v, w, sv-eps, err, sl_thresh)) / (err); // err = 2 * eps

  // double taa = t / (a*a);
  // double nnt = (1 + sv*sv * t);
  // double onnt = 1 / nnt;
  // double sqtonnt = sqrt(onnt);
  // double mexp = exp(0.5 * onnt * (sv*sv*a*a*w*w - 2*v*a*w - v*v*t));
  // double arg = a*a*w*w - 2*v*a*w*t + v*v*t*t;
  // double mult = onnt*onnt*onnt*onnt * sqtonnt * (nnt * (arg - t - 3*sv*sv*t*t) +
  //               sv*sv * (arg - t*nnt) * (arg - 5*t*nnt));
  // double sum_err;

  // if (taa > sl_thresh) { // use large-time
  //   mult /= a*a;
  //   sum_err = err / fabs(mult);
  //   if (sum_err < ERR_TOL_THRESH) sum_err = ERR_TOL_THRESH;
  //   int kl = kl_pdf(taa, sum_err);
  //   return mult * PI_CONST * large_sum(taa, w, kl);
  // } else { // use small-time
  //   mult *= SQRT_1_2PI * a / (t * sqrt(t));
  //   sum_err = err / fabs(mult);
  //   if (sum_err < ERR_TOL_THRESH) sum_err = ERR_TOL_THRESH;
  //   return mult * small_sum(taa, w, sum_err);
  // }
}
