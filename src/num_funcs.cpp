// Functions to determine the number of terms required for approximating the
// infinite sum in the density function

#include <Rcpp.h>
// [[Rcpp::plugins(cpp11)]]
#define _USE_MATH_DEFINES
#include <math.h>
#include "funcs.h"




//////////                                                           //////////
///////////////////////////////// Small Time //////////////////////////////////
//////////                                                           //////////


// BGK2014
int ks_BGK(double t, const double& a, const double& w, const double& eps)
{
  double u_eps, arg, k1;
  int k;
  t /= a * a; // rescale t for one-parameter form
  u_eps = min(-1.0, log(2 * M_PI * t * t * eps * eps)); // Safe bound for sqrt
  arg = -t * (u_eps - sqrt(-2 * u_eps - 2)); // sqrt(x) with x > 0
  k1 = (sqrt(2 * t) - w)/2;
  if (arg > 0) { // If arg > 0, set k2 and calculate k
    double k2 = (sqrt(arg) - w) / 2;
    k = ceil(std::max(k1, k2));
  }
  else { // Otherwise, we don't need k2
    k = ceil(k1);
  }
  return k;
}


// Navarro2009
int ks_Nav(double t, const double& a, const double& eps)
{
  t /= a * a; // rescale t for one-parameter form
  if (eps < 1 / (2 * sqrt(2 * M_PI * t))) { // if error threshold is set low enough
    double ks = 2 + sqrt(-2 * t * log(2 * eps * sqrt(2 * M_PI * t)));
    return ceil(std::max(ks, sqrt(t)+1)); // ensure boundary conditions are met
  }
  return 2; // else return minimal kappa for that case
}




//////////                                                           //////////
///////////////////////////////// Large Time //////////////////////////////////
//////////                                                           //////////


// Navarro2009
int kl_Nav(double t, const double& a, const double& eps)
{
  t /= a * a; // rescale t for one-parameter form
  if (eps < 1 / (M_PI * t)) { // error threshold is low enough
    double kl = sqrt(-2 * log(M_PI * t * eps) / (M_PI * M_PI *t));
    return ceil(std::max(kl, 1 / (M_PI * sqrt(t))));
  }
  return ceil(1 / (M_PI * sqrt(t))); // else set to boundary condition
}
