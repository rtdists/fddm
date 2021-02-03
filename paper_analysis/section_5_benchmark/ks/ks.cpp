// C++ functions for calculating the required number of terms in the truncated
// "small-time" sum. These functions are identical to the ones in
// `src/num_funcs.cpp`.

#include <Rcpp.h>
// [[Rcpp::plugins(cpp11)]]
using std::min;
using std::max;
using std::log;
using std::sqrt;
using std::ceil;



////////// Gondan et al 2014 ///////////////////////////////////////////////////
// [[Rcpp::export]]
int ks_Gon(const double& t, const double& w, const double& eps)
{
  float u_eps, arg, k1;
  int k;
  u_eps = min(-1.0, log(2 * M_PI * t*t * eps*eps)); // Safe bound for sqrt
  arg = -t * (u_eps - sqrt(-2 * u_eps - 2)); // sqrt(x) with x > 0
  k1 = (sqrt(2 * t) - w)/2;
  if (k1 > INT_MAX) return INT_MAX;
  if (arg > 0) { // If arg > 0, set k2 and calculate k
    float k2 = (sqrt(arg) - w) / 2;
    if (k2 > INT_MAX) return INT_MAX;
    k = ceil(max(k1, k2));
  }
  else { // Otherwise, we don't need k2
    k = ceil(k1);
  }
  return 2 * k + 1; // Convert *pairs* of terms to total *individual* terms
}


////////// Navarro and Fuss 2009 ///////////////////////////////////////////////
// [[Rcpp::export]]
int ks_Nav(const double& t, const double& w, const double& eps)
{ // note: w is not used
  if (eps < 1 / (2 * sqrt(2 * M_PI * t))) { // if error threshold is set low enough
    float ks = 2 + sqrt(-2 * t * log(2 * eps * sqrt(2 * M_PI * t)));
    float bc = sqrt(t) + 1; // boundary conditions
    if (ks > INT_MAX || bc > INT_MAX) return INT_MAX;
    return ceil(max(ks, bc)); // ensure boundary conditions are met
  }
  return 2; // else return minimal kappa for that case
}
