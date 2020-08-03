#include <Rcpp.h>
// [[Rcpp::plugins(cpp11)]]
#define _USE_MATH_DEFINES
#include <cmath>

using std::log;
using std::exp;
using std::sqrt;
using std::min;
using std::max;
using std::ceil;

static const double SQRT_2PI = sqrt(2 * M_PI);



// [[Rcpp::export]]
double ms(const double& t, const double& a, const double& v,
          const double& w, const double& sv)
{
  double mult_s;
  mult_s = a * exp((sv*sv * a*a * w*w - 2 * v * a * w - v*v * t)
                   / (2 + 2 * sv*sv * t))
             / ( (t * SQRT_2PI * sqrt(t)) * sqrt(1 + sv*sv * t) );
  return mult_s;
}


// [[Rcpp::export]]
double ml(const double& t, const double& a, const double& v,
          const double& w, const double& sv)
{
  // note: neither numm nor summ is used because there is only one option each
  double mult_l;
  mult_l = M_PI * exp((sv*sv * a*a * w*w - 2 * v * a * w - v*v * t)
                      / (2 + 2 * sv*sv * t))
                / (a*a * sqrt(1 + sv*sv * t));
  return mult_l;
}



// [[Rcpp::export]]
double ss17(const double& t, const double& a, const double& w,
            const int& ks)
{
  double gamma = -a*a / (2 * t);
  double rj;
  int j = 1;
  rj = j + 1 - w; // j is odd
  return rj * exp(gamma * rj*rj);
}

// [[Rcpp::export]]
double ss14(const double& t, const double& a, const double& w,
            const int& ks)
{
  // note: eps is not used
  double gamma = -a*a / (2 * t);
  int j = 1;
  return (2 * j + w) * exp(gamma * (2 * j + w) * (2 * j + w));
}


// [[Rcpp::export]]
double sl(const double& t, const double& a, const double& w,
          const int& kl)
{
  double gamma = -M_PI*M_PI * t / (2 * a*a);
  int j = 1;
  return j * sin(j * w * M_PI) * exp(gamma * j*j);
}


////////////////////////////////////////////////////////////////////////////////
// Num Funcs

// [[Rcpp::export]]
int ks_eps_17(const double& t, const double& a, const double& w, const double& eps)
{
  double gamma = -a*a / (2 * t);
  double sum = w * exp(gamma * w*w); // start at j=0 term
  double minterms = sqrt(t) / a; // minimum number of terms
  double rj = 1 - w;
  double term = rj * exp(gamma * rj*rj);
  int j = 0;
  while (j <= minterms) { // capture increasing terms
    j++;
    if (j % 2 == 1) { // j is odd
      rj = j + 1 - w;
      term = rj * exp(gamma * rj*rj);
      sum -= term;
    } else { // j is even
      rj = j + w;
      term = rj * exp(gamma * rj*rj);
      sum += term;
    }
  }
  while (fabs(term) > eps) { // at this point, odd (negative) term is greater
    j++;
    if (j % 2 == 1) { // j is odd
      rj = j + 1 - w;
      term = rj * exp(gamma * rj*rj);
      sum -= term;
    } else { // j is even
      rj = j + w;
      term = rj * exp(gamma * rj*rj);
      sum += term;
    }
  }
  return j + 1;
}

// [[Rcpp::export]]
int ks_eps_14(const double& t, const double& a, const double& w, const double& eps)
{
  double gamma = -a*a / (2 * t);
  double sum = w * exp(gamma * w*w); // start at j=0 term
  double minterms = sqrt(t) / a / 2; // minimum number of terms
  double pterm = 0;
  double nterm = sum;
  int j = 0;
  while (j <= minterms) { // capture increasing terms
    j++;
    nterm = (w - 2 * j) * exp(gamma * (w - 2 * j) * (w - 2 * j));
    pterm = (w + 2 * j) * exp(gamma * (w + 2 * j) * (w + 2 * j));
    sum += nterm + pterm;
  }
  while (fabs(pterm) > eps) { // at this point, the negative term is greater
    j++;
    nterm = (w - 2 * j) * exp(gamma * (w - 2 * j) * (w - 2 * j));
    pterm = (w + 2 * j) * exp(gamma * (w + 2 * j) * (w + 2 * j));
    sum += nterm + pterm;
  }
  return 2*j + 1;
}


// [[Rcpp::export]]
int ks_Gon(const double& t, const double& w, const double& eps)
{
  double u_eps, arg, k1;
  int k;
  u_eps = min(-1.0, log(2 * M_PI * t*t * eps*eps)); // Safe bound for sqrt
  arg = -t * (u_eps - sqrt(-2 * u_eps - 2)); // sqrt(x) with x > 0
  k1 = (sqrt(2 * t) - w)/2;
  if (arg > 0) { // If arg > 0, set k2 and calculate k
    double k2 = (sqrt(arg) - w) / 2;
    k = ceil(max(k1, k2));
  }
  else { // Otherwise, we don't need k2
    k = ceil(k1);
  }
  return 2*k + 1;
}


// [[Rcpp::export]]
int ks_Nav(const double& t, const double& w, const double& eps)
{
  // note: w is not used
  if (eps < 1 / (2 * sqrt(2 * M_PI * t))) { // if error threshold is set low enough
    double ks = 2 + sqrt(-2 * t * log(2 * eps * sqrt(2 * M_PI * t)));
    return ceil(max(ks, sqrt(t)+1)); // ensure boundary conditions are met
  }
  return 2; // else return minimal kappa for that case
}


// [[Rcpp::export]]
int kl_Nav(const double& t, const double& w, const double& eps)
{
  // note: w is not used
  if (eps < 1 / (M_PI * t)) { // error threshold is low enough
    double kl = sqrt(-2 * log(M_PI * t * eps) / (M_PI*M_PI * t));
    return ceil(max(kl, 1 / (M_PI * sqrt(t))));
  }
  return ceil(1 / (M_PI * sqrt(t))); // else set to boundary condition
}
