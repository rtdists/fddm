// Density functions for the Ratcliff diffusion model

#include <Rcpp.h>
// [[Rcpp::plugins(cpp11)]]
#define _USE_MATH_DEFINES
#include <math.h>




////////////////////////////////////////////////////////////////////////////////
///////////// Number of Terms Approximations

///// Small Time /////
// BGK2014
// [[Rcpp::export]]
int ks_BGK(double t, double a, double w, double eps) {
  double u_eps, arg, k1;
  int k;
  t /= a*a; // rescale t for one-parameter form
  u_eps = std::min(-1.0, log(2*M_PI*t*t*eps*eps)); // Safe bound for sqrt below
  arg = -t*(u_eps - sqrt(-2*u_eps - 2)); // sqrt(x) with x > 0
  k1 = (sqrt(2*t) - w)/2;
  if (arg > 0) { // If arg > 0, set k2 and calculate k
    double k2 = (sqrt(arg) - w)/2;
    k = ceil(std::max(k1, k2));
  }
  else { // Otherwise, we don't need k2
    k = ceil(k1);
  }
  return k;
}

// Navarro2009
// [[Rcpp::export]]
int ks_Nav(double t, double a, double eps)
{
  t /= a*a; // rescale t for one-parameter form
  if (eps < 1/(2 * sqrt(2*M_PI*t))) { // if error threshold is set low enough
    double ks = 2 + sqrt(-2*t*log(2*eps*sqrt(2*M_PI*t)));
    return ceil(std::max(ks, sqrt(t)+1)); // ensure boundary conditions are met
  }
  return 2; // else return minimal kappa for that case
}


///// Large Time /////
// Navarro2009
// [[Rcpp::export]]
int kl_Nav(double t, double a, double eps)
{
  t /= a*a; // rescale t for one-parameter form
  if (eps < 1/(M_PI*t)) { // error threshold is low enough
    double kl = sqrt(-2*log(M_PI*t*eps)/(M_PI*M_PI*t));
    return ceil(std::max(kl, 1/(M_PI*sqrt(t))));
  }
  return 1/(M_PI*sqrt(t)); // I don't think I need this actually
}




////////////////////////////////////////////////////////////////////////////////
///////////// Infinite Sum Approximations

///// Small Time /////
// BGK2014 style truncated sum
// [[Rcpp::export]]
double small_sum_2014(double t, double a, double w, int ks)
{
  double gamma = -a*a/(2*t);
  double f = w * exp(gamma*w*w); // start at j=0
  for (int j = ks; j > 0; j--) { // iterate through all ks
    f += (2*j + w) * exp(gamma*(2*j + w)*(2*j + w))
       - (2*j - w) * exp(gamma*(2*j - w)*(2*j - w));
  }
  return f;
}

// BGK2017 style truncated sum (uses ceil(k/2) as many terms as BGK2014)
// [[Rcpp::export]]
double small_sum_2017(double t, double a, double w, int ks)
{
  double gamma = -a*a/(2*t);
  double sum = w * exp(gamma*w*w); // start at j=0
  double rj;
  int j = 1;
  while (j <= 2*ks) { // start at j=1
    rj = j + 1 - w; // j is odd
    sum -= rj * exp(gamma*rj*rj);
    j++;
    rj = j + w; // j is even
    sum += rj * exp(gamma*rj*rj);
    j++;
  }
  return sum;
}

// term < eps BGK2014 style truncated sum, with minimum terms
// [[Rcpp::export]]
double small_sum_eps_14(double t, double a, double w, double mult, double eps)
{
  eps /= mult;
  double gamma = -a*a/(2*t);
  double sum = w * exp(gamma * w * w); // start at j=0 term
  double minterms = sqrt(t)/(2*a) - w/2; // minimum number of terms
  double pterm = 0;
  double nterm = sum;
  int j = 0;
  while (j < minterms) { // capture increasing terms
    j++;
    pterm = (w + 2*j) * exp(gamma * (w + 2*j) * (w + 2*j));
    nterm = (w - 2*j) * exp(gamma * (w - 2*j) * (w - 2*j));
    sum += pterm + nterm;
  }
  while (fabs(nterm) > eps) { // at this point, the negative term is greater
    j++;
    pterm = (w + 2*j) * exp(gamma * (w + 2*j) * (w + 2*j));
    nterm = (w - 2*j) * exp(gamma * (w - 2*j) * (w - 2*j));
    sum += pterm + nterm;
  }
  return sum;
}

// term < eps BGK2017 style truncated sum, with minimum terms
// [[Rcpp::export]]
double small_sum_eps_17(double t, double a, double w, double mult, double eps)
{
  eps /= mult;
  double gamma = -a*a/(2*t);
  double sum = w * exp(gamma * w * w); // start at j=0 term
  double minterms = sqrt(t)/(2*a) - w/2; // minimum number of terms
  double rj = 1 - w;
  double oterm = rj * exp(gamma*rj*rj);
  double eterm = 0;
  int j = 0;
  while (j < minterms) { // capture increasing terms
    j++;
    rj = j + 1 - w; // j is odd
    oterm = rj * exp(gamma*rj*rj);
    sum -= oterm;
    j++;
    rj = j + w; // j is even
    eterm = rj * exp(gamma*rj*rj);
    sum += eterm;
  }
  while (fabs(oterm) > eps) { // at this point, the odd (negative) term is greater
    j++;
    rj = j + 1 - w; // j is odd
    oterm = rj * exp(gamma*rj*rj);
    sum -= oterm;
    j++;
    rj = j + w; // j is even
    eterm = rj * exp(gamma*rj*rj);
    sum += eterm;
  }
  return sum;
}


///// Large Time /////
// Navarro2009 style truncated sum
// [[Rcpp::export]]
double large_sum_Nav(double t, double a, double w, int kl)
{
  double gamma = -M_PI*M_PI*t/(2*a*a);
  double sum = 0.0;
  for (int j = 1; j <= kl; j++) {
    sum += j * sin(j*w*M_PI) * exp(gamma*j*j);
  }
  return sum;
}




////////////////////////////////////////////////////////////////////////////////
///////////// Densities

///// Small Time /////
// Use term < eps BGK2014 style sum approximation, with minimum terms
// [[Rcpp::export]]
Rcpp::NumericVector fs_eps_2014(Rcpp::NumericVector rt, Rcpp::LogicalVector resp,
                                double v, double a, double w,
                                double t0=0.0, double eps=0.0)
{
  if (eps <= 0.0) { // set eps to default value (doesn't work in definition)
    eps = sqrt(DBL_EPSILON);
  }

  int n = rt.size(); // get number of values in rt
  Rcpp::NumericVector out(n);
  double t, mult;
  double vprime = -v;
  double wprime = 1 - w;

  for (int i = 0; i < n; i++) { // iterate through all values in rt and out
    t = rt[i] - t0; // subtract non-decisison time from the response time
    if (t <= 0) {
      out[i] = 0; // if t=0, return 0 instead of NaN
      continue;
    }

    if (resp[i] == 1) { // if response is "upper" then change parameter values
      mult = a * exp(-vprime*a*wprime - vprime*vprime*t/2) / sqrt(2*M_PI*t*t*t);
      out[i] = mult * small_sum_eps_14(t, a, wprime, mult, eps);
    } else { // else response is "lower"
      mult = a * exp(-v*a*w - v*v*t/2) / sqrt(2*M_PI*t*t*t);
      out[i] = mult * small_sum_eps_14(t, a, w, mult, eps);
    }
  }
  return out;
}

// Use term < eps BGK2017 style sum approximation, with minimum terms
// [[Rcpp::export]]
Rcpp::NumericVector fs_eps_2017(Rcpp::NumericVector rt, Rcpp::LogicalVector resp,
                                double v, double a, double w,
                                double t0=0.0, double eps=0.0)
{
  if (eps <= 0.0) { // set eps to default value (doesn't work in definition)
    eps = sqrt(DBL_EPSILON);
  }

  int n = rt.size(); // get number of values in rt
  Rcpp::NumericVector out(n);
  double t, mult;
  double vprime = -v;
  double wprime = 1 - w;

  for (int i = 0; i < n; i++) { // iterate through all values in rt and out
    t = rt[i] - t0; // subtract non-decisison time from the response time
    if (t <= 0) {
      out[i] = 0; // if t=0, return 0 instead of NaN
      continue;
    }

    if (resp[i] == 1) { // if response is "upper" then change parameter values
      mult = a * exp(-vprime*a*wprime - vprime*vprime*t/2) / sqrt(2*M_PI*t*t*t);
      out[i] = mult * small_sum_eps_17(t, a, wprime, mult, eps);
    } else { // else response is "lower"
      mult = a * exp(-v*a*w - v*v*t/2) / sqrt(2*M_PI*t*t*t);
      out[i] = mult * small_sum_eps_17(t, a, w, mult, eps);
    }
  }
  return out;
}

// Use Navarro2009 number of terms for 2014 style sum approximation
// [[Rcpp::export]]
Rcpp::NumericVector fs_Nav_2014(Rcpp::NumericVector rt, Rcpp::LogicalVector resp,
                                double v, double a, double w,
                                double t0=0.0, double eps=0.0)
{
  if (eps <= 0.0) { // set eps to default value (doesn't work in definition)
    eps = sqrt(DBL_EPSILON);
  }

  int n = rt.size(); // get number of values in rt
  Rcpp::NumericVector out(n);
  int ks;
  double t, mult;
  double vprime = -v;
  double wprime = 1 - w;

  for (int i = 0; i < n; i++) { // iterate through all values in rt and out
    t = rt[i] - t0; // subtract non-decisison time from the response time
    if (t <= 0) {
      out[i] = 0; // if t=0, return 0 instead of NaN
      continue;
    }

    ks = ks_Nav(t, a, eps); // get number of terms in sum approximation
    if (resp[i] == 1) { // if response is "upper" then change parameter values
      mult = a * exp(-vprime*a*wprime - vprime*vprime*t/2) / sqrt(2*M_PI*t*t*t);
      out[i] = mult * small_sum_2014(t, a, wprime, ks);
    } else { // else response is "lower"
      mult = a * exp(-v*a*w - v*v*t/2) / sqrt(2*M_PI*t*t*t);
      out[i] = mult * small_sum_2014(t, a, w, ks);
    }
  }
  return out;
}

// Use Navarro2009 number of terms for 2017 style sum approximation
// [[Rcpp::export]]
Rcpp::NumericVector fs_Nav_2017(Rcpp::NumericVector rt, Rcpp::LogicalVector resp,
                                double v, double a, double w,
                                double t0=0.0, double eps=0.0)
{
  if (eps <= 0.0) { // set eps to default value (doesn't work in definition)
    eps = sqrt(DBL_EPSILON);
  }

  int n = rt.size(); // get number of values in rt
  Rcpp::NumericVector out(n);
  int ks;
  double t, mult;
  double vprime = -v;
  double wprime = 1 - w;

  for (int i = 0; i < n; i++) { // iterate through all values in rt and out
    t = rt[i] - t0; // subtract non-decisison time from the response time
    if (t <= 0) {
      out[i] = 0; // if t=0, return 0 instead of NaN
      continue;
    }

    ks = ks_Nav(t, a, eps); // get number of terms in sum approximation
    if (resp[i] == 1) { // if response is "upper" then change parameter values
      mult = a * exp(-vprime*a*wprime - vprime*vprime*t/2) / sqrt(2*M_PI*t*t*t);
      out[i] = mult * small_sum_2017(t, a, wprime, ks);
    } else { // else response is "lower"
      mult = a * exp(-v*a*w - v*v*t/2) / sqrt(2*M_PI*t*t*t);
      out[i] = mult * small_sum_2017(t, a, w, ks);
    }
  }
  return out;
}

// Use BGK2014 number of terms for 2014 style sum approximation
// [[Rcpp::export]]
Rcpp::NumericVector fs_BGK_2014(Rcpp::NumericVector rt, Rcpp::LogicalVector resp,
                                double v, double a, double w,
                                double t0=0.0, double eps=0.0)
{
  if (eps <= 0.0) { // set eps to default value (doesn't work in definition)
    eps = sqrt(DBL_EPSILON);
  }

  int n = rt.size(); // get number of values in rt
  Rcpp::NumericVector out(n);
  int ks;
  double t, mult;
  double vprime = -v;
  double wprime = 1 - w;

  for (int i = 0; i < n; i++) { // iterate through all values in rt and out
    t = rt[i] - t0; // subtract non-decisison time from the response time
    if (t <= 0) {
      out[i] = 0; // if t=0, return 0 instead of NaN
      continue;
    }

    if (resp[i] == 1) { // if response is "upper" then change parameter values
      ks = ks_BGK(t, a, wprime, eps); // get number of terms in sum approximation
      mult = a * exp(-vprime*a*wprime - vprime*vprime*t/2) / sqrt(2*M_PI*t*t*t);
      out[i] = mult * small_sum_2014(t, a, wprime, ks);
    } else { // else response is "lower"
      ks = ks_BGK(t, a, w, eps); // get number of terms in sum approximation
      mult = a * exp(-v*a*w - v*v*t/2) / sqrt(2*M_PI*t*t*t);
      out[i] = mult * small_sum_2014(t, a, w, ks);
    }
  }
  return out;
}

// Use BGK2014 number of terms for 2017 style sum approximation
// [[Rcpp::export]]
Rcpp::NumericVector fs_BGK_2017(Rcpp::NumericVector rt, Rcpp::LogicalVector resp,
                                double v, double a, double w,
                                double t0=0.0, double eps=0.0)
{
  if (eps <= 0.0) { // set eps to default value (doesn't work in definition)
    eps = sqrt(DBL_EPSILON);
  }

  int n = rt.size(); // get number of values in rt
  Rcpp::NumericVector out(n);
  int ks;
  double t, mult;
  double vprime = -v;
  double wprime = 1 - w;

  for (int i = 0; i < n; i++) { // iterate through all values in rt and out
    t = rt[i] - t0; // subtract non-decisison time from the response time
    if (t <= 0) {
      out[i] = 0; // if t=0, return 0 instead of NaN
      continue;
    }

    if (resp[i] == 1) { // if response is "upper" then change parameter values
      ks = ks_BGK(t, a, wprime, eps); // get number of terms in sum approximation
      mult = a * exp(-vprime*a*wprime - vprime*vprime*t/2) / sqrt(2*M_PI*t*t*t);
      out[i] = mult * small_sum_2017(t, a, wprime, ks);
    } else { // else response is "lower"
      ks = ks_BGK(t, a, w, eps); // get number of terms in sum approximation
      mult = a * exp(-v*a*w - v*v*t/2) / sqrt(2*M_PI*t*t*t);
      out[i] = mult * small_sum_2017(t, a, w, ks);
    }
  }
  return out;
}

///// Large Time /////
// Use Navarro2009 number of terms for sum approximation
// [[Rcpp::export]]
Rcpp::NumericVector fl_Nav(Rcpp::NumericVector rt, Rcpp::LogicalVector resp,
                           double v, double a, double w,
                           double t0=0.0, double eps=0.0)
{
  if (eps <= 0.0) { // set eps to default value (doesn't work in definition)
    eps = sqrt(DBL_EPSILON);
  }

  int n = rt.size(); // get number of values in rt
  Rcpp::NumericVector out(n);
  double t, mult;
  int kl;
  double vprime = -v;
  double wprime = 1 - w;

  for (int i = 0; i < n; i++) { // iterate through all values in rt and out
    t = rt[i] - t0; // subtract non-decisison time from the response time
    if (t <= 0) {
      out[i] = 0; // if t=0, return 0 instead of NaN
      continue;
    }

    kl = kl_Nav(t=t, a=a, eps=eps);
    if (resp[i] == 1) { // if response is "upper" then change parameter values
      mult = M_PI * exp(-vprime*a*wprime - vprime*vprime*t/2) / (a*a);
      out[i] = mult * large_sum_Nav(t, a, wprime, kl);
    } else { // else response is "lower"
      mult = M_PI * exp(-v*a*w - v*v*t/2) / (a*a);
      out[i] = mult * large_sum_Nav(t, a, w, kl);
    }
  }
  return out;
}


///// Combined Small and Large Time /////
// ks = Navarro2009, 2014 style sum approximation, kl = Navarro2009
// [[Rcpp::export]]
Rcpp::NumericVector fb_Nav_Nav_14(Rcpp::NumericVector rt,
                                  Rcpp::LogicalVector resp,
                                  double v, double a, double w,
                                  double t0=0.0, double eps=0.0)
{
  if (eps <= 0.0) { // set eps to default value (doesn't work in definition)
    eps = sqrt(DBL_EPSILON);
  }

  int n = rt.size(); // get number of values in rt
  Rcpp::NumericVector out(n);
  double t, mult;
  int ks, kl;
  double vprime = -v;
  double wprime = 1 - w;

  for (int i = 0; i < n; i++) { // iterate through all values in rt and out
    t = rt[i] - t0; // subtract non-decisison time from the response time
    if (t <= 0) {
      out[i] = 0; // if t=0, return 0 instead of NaN
      continue;
    }

    ks = ks_Nav(t, a, eps);
    kl = kl_Nav(t, a, eps);
    if (resp[i] == 1) { // if response is "upper" then change parameter values
      if(ks < kl) {
        mult = a * exp(-vprime*a*wprime - vprime*vprime*t/2) / sqrt(2*M_PI*t*t*t);
        out[i] = mult * small_sum_2014(t, a, wprime, ks);
      }
      else {
        mult = M_PI * exp(-vprime*a*wprime - vprime*vprime*t/2) / (a*a);
        out[i] = mult * large_sum_Nav(t, a, wprime, kl);
      }
    } else { // else response is "lower"
      if(ks < kl) {
        mult = a * exp(-v*a*w - v*v*t/2) / sqrt(2*M_PI*t*t*t);
        out[i] = mult * small_sum_2014(t, a, w, ks);
      }
      else {
        mult = M_PI * exp(-v*a*w - v*v*t/2) / (a*a);
        out[i] = mult * large_sum_Nav(t, a, w, kl);
      }
    }
  }
  return out;
}

// ks = Navarro2009, 2017 style sum approximation, kl = Navarro2009
// [[Rcpp::export]]
Rcpp::NumericVector fb_Nav_Nav_17(Rcpp::NumericVector rt,
                                  Rcpp::LogicalVector resp,
                                  double v, double a, double w,
                                  double t0=0.0, double eps=0.0)
{
  if (eps <= 0.0) { // set eps to default value (doesn't work in definition)
    eps = sqrt(DBL_EPSILON);
  }

  int n = rt.size(); // get number of values in rt
  Rcpp::NumericVector out(n);
  double t, mult;
  int ks, kl;
  double vprime = -v;
  double wprime = 1 - w;

  for (int i = 0; i < n; i++) { // iterate through all values in rt and out
    t = rt[i] - t0; // subtract non-decisison time from the response time
    if (t <= 0) {
      out[i] = 0; // if t=0, return 0 instead of NaN
      continue;
    }

    ks = ks_Nav(t, a, eps);
    kl = kl_Nav(t, a, eps);
    if (resp[i] == 1) { // if response is "upper" then change parameter values
      if(ks < kl) {
        mult = a * exp(-vprime*a*wprime - vprime*vprime*t/2) / sqrt(2*M_PI*t*t*t);
        out[i] = mult * small_sum_2017(t=t, a=a, w=wprime, ks=ks);
      }
      else {
        mult = M_PI * exp(-vprime*a*wprime - vprime*vprime*t/2) / (a*a);
        out[i] = mult * large_sum_Nav(t, a, wprime, kl);
      }
    } else { // else response is "lower"
      if(ks < kl) {
        mult = a * exp(-v*a*w - v*v*t/2) / sqrt(2*M_PI*t*t*t);
        out[i] = mult * small_sum_2017(t, a, w, ks);
      }
      else {
        mult = M_PI * exp(-v*a*w - v*v*t/2) / (a*a);
        out[i] = mult * large_sum_Nav(t, a, w, kl);
      }
    }
  }
  return out;
}

// ks = BGK2014, 2014 style sum approximation, kl = Navarro2009
// [[Rcpp::export]]
Rcpp::NumericVector fb_BGK_Nav_14(Rcpp::NumericVector rt,
                                  Rcpp::LogicalVector resp,
                                  double v, double a, double w,
                                  double t0=0.0, double eps=0.0)
{
  if (eps <= 0.0) { // set eps to default value (doesn't work in definition)
    eps = sqrt(DBL_EPSILON);
  }

  int n = rt.size(); // get number of values in rt
  Rcpp::NumericVector out(n);
  double t, mult;
  int ks, kl;
  double vprime = -v;
  double wprime = 1 - w;

  for (int i = 0; i < n; i++) { // iterate through all values in rt and out
    t = rt[i] - t0; // subtract non-decisison time from the response time
    if (t <= 0) {
      out[i] = 0; // if t=0, return 0 instead of NaN
      continue;
    }

    kl = kl_Nav(t, a, eps);
    if (resp[i] == 1) { // if response is "upper" then change parameter values
      ks = ks_BGK(t, a, wprime, eps);
      if(ks < kl) {
        mult = a * exp(-vprime*a*wprime - vprime*vprime*t/2) / sqrt(2*M_PI*t*t*t);
        out[i] = mult * small_sum_2014(t, a, wprime, ks);
      }
      else {
        mult = M_PI * exp(-vprime*a*wprime - vprime*vprime*t/2) / (a*a);
        out[i] = mult * large_sum_Nav(t, a, wprime, kl);
      }
    } else { // else response is "lower"
      ks = ks_BGK(t, a, w, eps);
      if(ks < kl) {
        mult = a * exp(-v*a*w - v*v*t/2) / sqrt(2*M_PI*t*t*t);
        out[i] = mult * small_sum_2014(t, a, w, ks);
      }
      else {
        mult = M_PI * exp(-v*a*w - v*v*t/2) / (a*a);
        out[i] = mult * large_sum_Nav(t, a, w, kl);
      }
    }
  }
  return out;
}

// ks = BGK2014, 2017 style sum approximation, kl = Navarro2009
// [[Rcpp::export]]
Rcpp::NumericVector fb_BGK_Nav_17(Rcpp::NumericVector rt,
                                  Rcpp::LogicalVector resp,
                                  double v, double a, double w,
                                  double t0=0.0, double eps=0.0)
{
  if (eps <= 0.0) { // set eps to default value (doesn't work in definition)
    eps = sqrt(DBL_EPSILON);
  }

  int n = rt.size(); // get number of values in rt
  Rcpp::NumericVector out(n);
  double t, mult;
  int ks, kl;
  double vprime = -v;
  double wprime = 1 - w;

  for (int i = 0; i < n; i++) { // iterate through all values in rt and out
    t = rt[i] - t0; // subtract non-decisison time from the response time
    if (t <= 0) {
      out[i] = 0; // if t=0, return 0 instead of NaN
      continue;
    }

    kl = kl_Nav(t, a, eps);
    if (resp[i] == 1) { // if response is "upper" then change parameter values
      ks = ks_BGK(t, a, wprime, eps);
      if(ks < kl) {
        mult = a * exp(-vprime*a*wprime - vprime*vprime*t/2) / sqrt(2*M_PI*t*t*t);
        out[i] = mult * small_sum_2017(t, a, wprime, ks);
      }
      else {
        mult = M_PI * exp(-vprime*a*wprime - vprime*vprime*t/2) / (a*a);
        out[i] = mult * large_sum_Nav(t, a, wprime, kl);
      }
    } else { // else response is "lower"
      ks = ks_BGK(t, a, w, eps);
      if(ks < kl) {
        mult = a * exp(-v*a*w - v*v*t/2) / sqrt(2*M_PI*t*t*t);
        out[i] = mult * small_sum_2017(t, a, w, ks);
      }
      else {
        mult = M_PI * exp(-v*a*w - v*v*t/2) / (a*a);
        out[i] = mult * large_sum_Nav(t, a, w, kl);
      }
    }
  }
  return out;
}


///// Variable Drift Rate (Small Time) /////
// [[Rcpp::export]]
Rcpp::NumericVector fs_vary(Rcpp::NumericVector rt, Rcpp::LogicalVector resp,
                            double v, double a, double w, double sv=0.0,
                            double t0=0.0, double eps=0.0)
{
  if (eps <= 0.0) { // set eps to default value (doesn't work in definition)
    eps = sqrt(DBL_EPSILON);
  }

  int n = rt.size(); // get number of values in rt
  Rcpp::NumericVector out(n);
  double t, mult;
  double vprime = -v;
  double wprime = 1 - w;

  if (sv < 0.05) { // set sv=0 and use constant drift rate method
    for (int i = 0; i < n; i++) { // iterate through all values in rt and out
      t = rt[i] - t0; // subtract non-decisison time from the response time
      if (t <= 0) {
        out[i] = 0; // if t=0, return 0 instead of NaN
        continue;
      }

      if (resp[i] == 1) { // if response is "upper" then change parameter values
        mult = a * exp(-vprime*a*wprime - vprime*vprime*t/2) / sqrt(2*M_PI*t*t*t);
        out[i] = mult * small_sum_eps_17(t, a, wprime, mult, eps);
      } else { // else response is "lower"
        mult = a * exp(-v*a*w - v*v*t/2) / sqrt(2*M_PI*t*t*t);
        out[i] = mult * small_sum_eps_17(t, a, w, mult, eps);
      }
    }
  }
  else { // use variable drift rate (changes mult)
    for (int i = 0; i < n; i++) { // iterate through all values in rt and out
      t = rt[i] - t0; // subtract non-decisison time from the response time
      if (t <= 0) {
        out[i] = 0; // if t=0, return 0 instead of NaN
        continue;
      }

      if (resp[i] == 1) { // if response is "upper" then change parameter values
        mult = a * exp((-vprime*vprime*t - 2*vprime*a*wprime + sv*a*a*wprime*wprime)
               / (2 + 2*sv*t)) / (sqrt(t*t*t + sv*t*t*t*t));
        out[i] = mult * small_sum_eps_17(t, a, wprime, mult, eps);
      } else { // else response is "lower"
        mult = a * exp((-v*v*t - 2*v*a*w + sv*a*a*w*w)
               / (2 + 2*sv*t)) / (sqrt(t*t*t + sv*t*t*t*t));
        out[i] = mult * small_sum_eps_17(t, a, w, mult, eps);
      }
    }
  }

  return out;
}
