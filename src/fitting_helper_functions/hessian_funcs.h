// Functions to evaluate the gradient of the DDM PDF
// included from src/fitting_helper_functions/class_methods.h



static const double ERR_TOL_THRESH = 1e-300; // near minimum value of a double
static const double SV_THRESH = 0; // threshold for using variable drift rate

static const double FOUR_NINTHS = 0.444444444444444444444;
static const double SQRT_2 = 1.41421356237309504880;
static const double SQRT_3 = sqrt(3);
static const double SQRT_5 = sqrt(5);
static const double LOG_2 = 0.693147180559945309417;
static const double LOG_4RT2 = log(4.0 * SQRT_2);
static const double LOG_5_112 = log(5.0/112.0);

static const double PI_CONST = 3.14159265358979323846; // define pi like C++
static const double LOG_PI = log(PI_CONST);
static const double LOG_2PI_2 = 0.5 * log(2 * PI_CONST);
static const double O_PI = 0.318309886183790671538;
static const double PI2 = PI_CONST * PI_CONST;
static const double PI3 = PI_CONST * PI_CONST * PI_CONST;
static const double SQRT_2PI = sqrt(2 * PI_CONST);
static const double SQRT_1_2PI = 1 / SQRT_2PI;
static const double SQRT_2_PI = sqrt(2 * O_PI);
static const double SQRT_2_1_PI = SQRT_2 * O_PI;

#include "num_funcs.h"
#include "sum_funcs.h"
#include "gradient_funcs.h"



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
{ // only uses small-time
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
    m2 = mexp * sqtonnt * PI3*PI2 * taa*taa / (a*a*a*a);
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
{
  Rcpp::Function dt_grad("numDeriv::grad");
  return dt_grad(_[func]=);
}

double dt02(const double& t, const double& v, const double& a, const double& w,
            const double& sv, const double& err, const double& sl_thresh)
{ // only uses small-time
  double taa = t / (a*a);
  double ot = 1 / t;
  double sqtot = sqrt(ot);
  double nnt = 1 + sv*sv * t;
  double nnt34 = 3 + 4 * sv*sv * t;
  double nnt38 = 3 + 8 * sv*sv * t;
  double nnt58 = 5 + 8 * sv*sv * t;
  double onnt = 1 / nnt;
  double sqtonnt = sqrt(onnt);
  double nawv = sv*sv*a*w - v;
  double nawvaw = sv*sv*a*a*w*w - v*a*w;
  double naw2vaw = sv*sv*a*a*w*w - 2*v*a*w;
  double arg = naw2vaw - v*v*t;
  double mexp = exp(0.5 * onnt * arg);
  double m0, m1, m2, sum_err0, sum_err1, sum_err2;
  
  // some decision between small-time and large-time
  if (taa > sl_thresh) { // use large-time
    m0 = 0.25 * mexp * onnt*onnt*onnt*onnt*sqtonnt / (a*a) * 
         (v*v*v*v + sv*sv*v*v*nnt + 3*sv*sv*sv*sv*nnt +
          sv*sv*naw2vaw*(2*v*v + 6*sv*sv*nnt + sv*sv*naw2vaw));
    m1 = -mexp * onnt*onnt*sqtonnt * (sv*sv*naw2vaw + (sv*sv + v*v)*nnt) /
         (a*a);
    m2 = 0.25 * mexp * PI3*PI2 * sqtonnt / (a*a*a*a*a*a);
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
    // m0 = -0.25 * mexp * SQRT_1_2PI * a * ot*ot*ot*sqtot *
    //      onnt*onnt*onnt*onnt*sqtonnt *
    //      (t*nnt*nnt *
    //         (14*sv*sv - 3*v*v - 3*sv*sv*nawvaw - nnt34*(5*sv*sv*sv*sv + v*v)) -
    //       t*t*nnt *
    //         (5*sv*sv*v*v + v*v*v*v + sv*sv*(5*sv*sv*nawvaw + v*v*naw2vaw)) -
    //       sv*sv*t*arg * (nnt*nnt34 + v*v*t + sv*sv*t*naw2vaw) +
    //       4*sv*sv*sv*sv*t*t*nnt*nnt - 5*nnt*nnt*nnt*nnt34);
    m0 = -0.25 * mexp * SQRT_1_2PI * a * ot*ot*ot*sqtot *
         onnt*onnt*onnt*onnt*sqtonnt *
         (nnt*( sv*sv*t*(8*nnt*nnt - 8*v*v - naw2vaw*nnt38) - 3*v*v*t - nnt*nnt34*nnt58) +
          t*nawv*nawv*(nnt*nnt34 + v*v*t + sv*sv*t*naw2vaw)
         );
    m1 = -0.25 * mexp * SQRT_1_2PI * a*a*a * ot*ot*ot*ot*sqtot *
         onnt*onnt*sqtonnt *
         (nnt*nnt34 + v*v*t + sv*sv*t*naw2vaw + ot*ot*ot*sqtot * sqtonnt *
          (nnt*(14 + sv*sv*t + v*t) + sv*sv*t * arg));
    m2 = 0.125 * mexp * SQRT_1_2PI * a*a*a*a*a*a*a *
         ot*ot*ot*ot*ot*ot*ot*sqtot * sqtonnt;
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
{ // only uses small-time
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


// [[Rcpp::export]]
vector<double> wrap(const vector<double>& t, const vector<double>& v,
                    const vector<double>& a, const vector<double>& w,
                    const vector<double>& sv, const double& err,
                    const double& sl_thresh, const int& par)
{
  int n = t.size();
  vector<double> out(n);
  if (par == 1) {
    for (int i = 0; i < n; i++) {
      out[i] = dv2(t[i], v[i], a[i], w[i], sv[i], err, sl_thresh);
    }
  } else if (par == 2) {
    for (int i = 0; i < n; i++) {
      out[i] = da2(t[i], v[i], a[i], w[i], sv[i], err, sl_thresh);
    }
  } else if (par == 3) {
    for (int i = 0; i < n; i++) {
      out[i] = dt02(t[i], v[i], a[i], w[i], sv[i], err, sl_thresh);
    }
  } else if (par == 4) {
    for (int i = 0; i < n; i++) {
      out[i] = dw2(t[i], v[i], a[i], w[i], sv[i], err, sl_thresh);
    }
  } else if (par == 5) {
    for (int i = 0; i < n; i++) {
      out[i] = dsv2(t[i], v[i], a[i], w[i], sv[i], err, sl_thresh);
    }
  }
  return out;
}

// [[Rcpp::export]]
vector<double> wrap1(const vector<double>& t, const vector<double>& v,
                    const vector<double>& a, const vector<double>& w,
                    const vector<double>& sv, const double& err,
                    const double& sl_thresh, const int& par)
{
  int n = t.size();
  vector<double> out(n);
  if (par == 1) {
    for (int i = 0; i < n; i++) {
      out[i] = dv(t[i], v[i], a[i], w[i], sv[i], err, sl_thresh);
    }
  } else if (par == 2) {
    for (int i = 0; i < n; i++) {
      out[i] = da(t[i], v[i], a[i], w[i], sv[i], err, sl_thresh);
    }
  } else if (par == 3) {
    for (int i = 0; i < n; i++) {
      out[i] = dt0(t[i], v[i], a[i], w[i], sv[i], err, sl_thresh);
    }
  } else if (par == 4) {
    for (int i = 0; i < n; i++) {
      out[i] = dw(t[i], v[i], a[i], w[i], sv[i], err, sl_thresh);
    }
  } else if (par == 5) {
    for (int i = 0; i < n; i++) {
      out[i] = dsv(t[i], v[i], a[i], w[i], sv[i], err, sl_thresh);
    }
  }
  return out;
}
