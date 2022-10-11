// Functions to evaluate the partial derivatives of the DDM PDF
// (for each parameter and the response time: t, a, v, w, sv)



//----------------- Regular (non-log) ----------------------------------------//

double pdf_dt(const double& t, const double& resp, const double& v,
              const double& a, const double& w, const double& sv,
              const double& err, const double& sl_thresh)
{ // note: resp is not used
  double taa = t / (a*a);
  double nnt = 1 + sv*sv * t;
  double sqtnnt = sqrt(nnt);
  double nawvaw = sv*sv * a*a * w*w - 2 * v * a * w;
  double mexp = exp(0.5 * (nawvaw - v*v * t) / nnt);

  if (taa > sl_thresh) { // use large-time
    double m1 = -0.5 * mexp * (v*v + sv*sv * (nawvaw + nnt)) /
                (a*a * nnt*nnt * sqtnnt);
    double m2 = mexp / (a*a * sqtnnt);
    int kl1 = kl_Nav(t / (a*a), 0.5 * err / fabs(m1));
    int kl2 = kl_Har(t / (a*a), t, 0.5 * err / fabs(m2));
    return m1 * PI_CONST * sum_large(taa, w, kl1) -
           0.5 * m2 * PI_CONST*PI_CONST*PI_CONST / (a*a) *
           sum_large_d(taa, w, kl2);
  } else { // use small-time
    double m1 = -0.5 * mexp * SQRT_1_2PI * a *
                (nnt * (3 + 4 * sv*sv * t) + v*v * t + sv*sv*t * nawvaw) /
                (t*t * sqrt(t) * nnt*nnt * sqtnnt);
    double m2 = 0.5 * mexp * SQRT_1_2PI * a*a*a / (t*t*t * sqrt(t) * sqtnnt);
    return m1 * sum_small(taa, w, 0.5 * err / fabs(m1)) +
           m2 * sum_small_d(taa, w, 0.5 * err / fabs(m2));
  }
}



double pdf_dt0(const double& t, const double& resp, const double& v,
               const double& a, const double& w, const double& sv,
               const double& err, const double& sl_thresh)
{ // note: resp is not used; also, this is just the negative of dt
  double taa = t / (a*a);
  double nnt = 1 + sv*sv * t;
  double sqtnnt = sqrt(nnt);
  double nawvaw = sv*sv * a*a * w*w - 2 * v * a * w;
  double mexp = exp(0.5 * (nawvaw - v*v * t) / nnt);

  if (taa > sl_thresh) { // use large-time
    double m1 = 0.5 * mexp * (v*v + sv*sv * (nawvaw + nnt)) /
                (a*a * nnt*nnt * sqtnnt);
    double m2 = -mexp / (a*a * sqtnnt);
    int kl1 = kl_Nav(t / (a*a), 0.5 * err / fabs(m1));
    int kl2 = kl_Har(t / (a*a), t, 0.5 * err / fabs(m2));
    return m1 * PI_CONST * sum_large(taa, w, kl1) -
           0.5 * m2 * PI_CONST*PI_CONST*PI_CONST / (a*a) *
           sum_large_d(taa, w, kl2);
  } else { // use small-time
    double m1 = 0.5 * mexp * SQRT_1_2PI * a *
                (nnt * (3 + 4 * sv*sv * t) + v*v * t + sv*sv*t * nawvaw) /
                (t*t * sqrt(t) * nnt*nnt * sqtnnt);
    double m2 = -0.5 * mexp * SQRT_1_2PI * a*a*a / (t*t*t * sqrt(t) * sqtnnt);
    return m1 * sum_small(taa, w, 0.5 * err / fabs(m1)) +
           m2 * sum_small_d(taa, w, 0.5 * err / fabs(m2));
  }
}



double pdf_dv(const double& t, const double& resp, const double& v,
              const double& a, const double& w, const double& sv,
              const double& err, const double& sl_thresh)
{ // note: resp = {1: "lower", 2: "upper"}
  double taa = t / (a*a);
  double nnt = 1 / (1 + sv*sv * t);
  double sqtnnt = sqrt(nnt);
  double arg = a * w + v * t;
  double mexp = exp(0.5 * (sv*sv * a*a * w*w - 2 * v * a * w - v*v * t) * nnt);
  double mult = (resp < 1.5) ? 1 : -1; // chain rule for "upper" response

  if (taa > sl_thresh) { // use large-time
    mult *= -mexp * arg * nnt * sqtnnt / (a*a);
    int kl = kl_Nav(taa, err / fabs(mult));
    return mult * PI_CONST * sum_large(taa, w, kl);
  } else { // use small-time
    mult *= -mexp * arg * a * SQRT_1_2PI * nnt * sqtnnt /
            (t * sqrt(t));
    return mult * sum_small(taa, w, err / fabs(mult));
  }
}



double pdf_da(const double& t, const double& resp, const double& v,
              const double& a, const double& w, const double& sv,
              const double& err, const double& sl_thresh)
{ // note: resp is not used
  double taa = t / (a*a);
  double nnt = 1 + sv*sv * t;
  double sqtnnt = sqrt(nnt);
  double naw = sv*sv * a*a * w*w;
  double vaw = v * a * w;
  double vvt = v*v * t;
  double mexp = exp(0.5 * (naw - 2 * vaw - vvt) / nnt);

  if (taa > sl_thresh) { // use large-time
    double m1 = mexp * (naw - vaw - 2 * nnt) / (a*a*a * nnt * sqtnnt);
    double m2 = mexp / (a*a * sqtnnt);
    int kl1 = kl_Nav(t / (a*a), 0.5 * err / fabs(m1));
    int kl2 = kl_Har(t / (a*a), t, 0.5 * err / fabs(m2));
    return m1 * PI_CONST * sum_large(taa, w, kl1) +
           m2 * PI_CONST*PI_CONST*PI_CONST * t / (a*a*a) *
           sum_large_d(taa, w, kl2);
  } else { // use small-time
    double m1 = mexp * (naw - vaw + nnt) * SQRT_1_2PI /
                (t * sqrt(t) * nnt * sqtnnt);
    double m2 = -mexp * a*a * SQRT_1_2PI / (t*t * sqrt(t) * sqtnnt);
    return m1 * sum_small(taa, w, 0.5 * err / fabs(m1)) +
           m2 * sum_small_d(taa, w, 0.5 * err / fabs(m2));
  }
}



double pdf_dw(const double& t, const double& resp, const double& v,
              const double& a, const double& w, const double& sv,
              const double& err, const double& sl_thresh)
{ // note: resp = {1: "lower", 2: "upper"}
  double out;
  double taa = t / (a*a);
  double onnt = 1 / (1 + sv*sv * t);
  double osqtnnt = sqrt(onnt);
  double arg = sv*sv * a * w - v;
  double mexp = exp(0.5 * onnt * (sv*sv * a*a * w*w - 2 * v * a * w - v*v * t));

  // first part
  double m1 = mexp * onnt * osqtnnt * arg / a;
  int kl1 = kl_Nav(t / (a*a), 0.5 * err / fabs(m1));
  if (kl1 <= sl_thresh) { // use large-time (same heuristic from dfddm)
    out = m1 * PI_CONST * sum_large(taa, w, kl1);
  } else { // use small-time
    m1 = mexp * arg * a*a * SQRT_1_2PI * onnt * osqtnnt / (t * sqrt(t));
    out = m1 * sum_small(taa, w, 0.5 * err / fabs(m1));
  }

  // second part
  double m2 = mexp * osqtnnt / (a*a);
  int kl2 = kl_Har_w(t / (a*a), t, 0.5 * err / fabs(m2));
  int ks2 = ks_Har_w(t / (a*a), w, 0.5 * err / fabs(m2));
  if (kl2 < 2 * ks2) { // use large-time
    out += m2 * PI_CONST*PI_CONST * sum_large_d_w(taa, w, kl2);
  } else { // use small time
    out += m2 * a*a*a * SQRT_1_2PI / (t * sqrt(t)) * sum_small_d_w(taa, w, ks2);
  }

  return (resp < 1.5) ? out : -out; // chain rule for "upper" response
}



double pdf_dsv(const double& t, const double& resp, const double& v,
               const double& a, const double& w, const double& sv,
               const double& err, const double& sl_thresh)
{ // note: resp is not used
  if (sv <= 0) { // check that the derivative actually makes sense
    warning("dsv_dfddm warning: function parameter 'sv' = 0.0; the derivative does not make sense; returning a NaN.");
    return std::numeric_limits<double>::quiet_NaN(); // return a NaN
  }
  double taa = t / (a*a);
  double nnt = 1 + sv*sv * t;
  double sqtnnt = sqrt(nnt);
  double arg = a*a * w*w + 2 * v * a * w * t + v*v * t*t - t * nnt;
  double mexp = exp(0.5 * (sv*sv * a*a * w*w - 2 * v * a * w - v*v * t) / nnt);
  double mult;

  if (taa > sl_thresh) { // use large-time
    mult = mexp * arg * sv / (a*a * nnt*nnt * sqtnnt);
    int kl = kl_Nav(taa, err / fabs(mult));
    return mult * PI_CONST * sum_large(taa, w, kl);
  } else { // use small-time
    mult = mexp * arg * sv * a * SQRT_1_2PI / (t * sqrt(t) * nnt*nnt * sqtnnt);
    return mult * sum_small(taa, w, err / fabs(mult));
  }
}



//----------------- Logged Versions ------------------------------------------//

// double pdf_dt_log(const double& t, const double& resp, const double& a,
//                   const double& v, const double& w, const double& sv,
//                   const double& err, const double& sl_thresh)
// {
//   double taa = t / (a*a);
//   double nnt = 1 + sv*sv * t;
//   double nawvaw = sv*sv * a*a * w*w - 2 * v * a * w;

//   if (taa > sl_thresh) { // use large-time
//     double m = 0.5 * PI2 / a;
//     int kl3 = kl_Har(taa, t, err / m);
//     int kl1 = kl_Nav(taa, err * m);
//     return -m * sum_large_d(taa, w, kl3) / sum_large(taa, w, kl1) -
//            0.5 * (sv*sv * nnt + v*v + sv*sv * nawvaw) / (nnt);
//   } else { // use small-time
//     double m = 0.5 * a / taa;
//     return m * sum_small_d(taa, w, err / m) / sum_small(taa, w, err * m) -
//            0.5 * (3 * nnt*nnt + sv*sv * t * nnt + v * t + sv*sv * t * nawvaw) /
//            (t * nnt*nnt);
//   }
// }



// double pdf_dt0_log(const double& t, const double& resp, const double& a,
//                    const double& v, const double& w, const double& sv,
//                    const double& err, const double& sl_thresh)
// { // note: resp is not used
//   double taa = t / (a*a);
//   double nnt = 1 + sv*sv * t;
//   double nawvaw = sv*sv * a*a * w*w - 2 * v * a * w;

//   if (taa > sl_thresh) { // use large-time
//     double m = 0.5 * PI2 / a;
//     int kl3 = kl_Har(taa, t, err / m);
//     int kl1 = kl_Nav(taa, err * m);
//     return m * sum_large_d(taa, w, kl3) / sum_large(taa, w, kl1) +
//            0.5 * (sv*sv * nnt + v*v + sv*sv * nawvaw) / (nnt);
//   } else { // use small-time
//     double m = 0.5 * a / taa;
//     return -m * sum_small_d(taa, w, err / m) / sum_small(taa, w, err * m) +
//            0.5 * (3 * nnt*nnt + sv*sv * t * nnt + v * t + sv*sv * t * nawvaw) /
//            (t * nnt*nnt);
//   }
// }



// double pdf_da_log(const double& t, const double& resp, const double& a,
//                   const double& v, const double& w, const double& sv,
//                   const double& err, const double& sl_thresh)
// {
//   double taa = t / (a*a);
//   double nnt = 1 + sv*sv * t;
//   double nawvaw = sv*sv * a*a * w*w - v * a * w;

//   if (taa > sl_thresh) { // use large-time
//     double m = PI2 * taa / a;
//     int kl3 = kl_Har(taa, t, err / m);
//     int kl1 = kl_Nav(taa, err * m);
//     return -m * sum_large_d(taa, w, kl3) / sum_large(taa, w, kl1) +
//            (nawvaw - 2 * nnt) / (a * nnt);
//   } else { // use small-time
//     double m = a / t;
//     return -m * sum_small_d(taa, w, err / m) / sum_small(taa, w, err * m) +
//            (nawvaw + nnt) / (a * nnt);
//   }
// }



// double pdf_dv_log(const double& t, const double& resp, const double& a,
//                   const double& v, const double& w, const double& sv,
//                   const double& err, const double& sl_thresh)
// {
//   return -1 * (a * w + v * t) / (1 + sv*sv * t);
// }



// double pdf_dw_log(const double& t, const double& resp, const double& a,
//                   const double& v, const double& w, const double& sv,
//                   const double& err, const double& sl_thresh)
// {
//   double taa = t / (a*a);
//   double nnt = 1 + sv*sv * t;
//   double arg = (sv*sv * a*a * w - v * a) / nnt;

//   int kl3 = kl_Har_w(t / (a*a), t, err / PI_CONST);
//   int kl1 = kl_Nav(t / (a*a), err * PI_CONST);
//   int ks3 = ks_Har_w(t / (a*a), w, err);

//   if (kl3 + kl1 < ks3 + 1) { // use large-time (same heuristic from dfddm)
//     return PI_CONST * sum_large_d_w(taa, w, kl3) / sum_large(taa, w, kl1) + arg;
//   } else { // use small-time
//     return sum_small_d_w(taa, w, ks3) / sum_small(taa, w, err) + arg;
//   }
// }



// double pdf_dsv_log(const double& t, const double& resp, const double& a,
//                    const double& v, const double& w, const double& sv,
//                    const double& err, const double& sl_thresh)
// {
//   if (sv <= 0) { // check that the derivative actually makes sense
//     warning("dsv_dfddm warning: function parameter 'sv' = 0.0; the derivative does not make sense; returning a NaN.");
//     return std::numeric_limits<double>::quiet_NaN(); // return a NaN
//   }
//   double nnt = 1 + sv*sv * t;
//   return -sv * t * (sv*sv * a*a * w*w - 2 * v * a * w - v*v * t + nnt) /
//          (nnt*nnt);
// }
