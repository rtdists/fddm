// C++ functions for calculating the required number of terms in the truncated
// "small-time" sum. These functions are based on the ones in
// `src/num_funcs.cpp`, but they accept a vector for the `t` input and loop over
// it. This is done to add more computation time in C++ so that the
// `microbenchmark` results contain proportionally less overhead from R.

#include <Rcpp.h>
// [[Rcpp::plugins(cpp11)]]
using std::min;
using std::max;
using std::log;
using std::sqrt;
using std::ceil;
using Rcpp::IntegerVector;
using Rcpp::NumericVector;



////////// Gondan et al 2014 ///////////////////////////////////////////////////
// [[Rcpp::export]]
IntegerVector ks_Gon(const NumericVector& t, const double& w, const double& eps)
{
  int n = t.length();
  IntegerVector out(n);
  double tt;
  float u_eps, arg, k1;
  int k;

  for (int i = 0; i < n; i++) {
    tt = t[i];

    u_eps = min(-1.0, log(2 * M_PI * tt*tt * eps*eps)); // Safe bound for sqrt
    arg = -tt * (u_eps - sqrt(-2 * u_eps - 2)); // sqrt(x) with x > 0
    k1 = (sqrt(2 * tt) - w)/2;
    if (k1 > INT_MAX) return INT_MAX;
    if (arg > 0) { // If arg > 0, set k2 and calculate k
      float k2 = (sqrt(arg) - w) / 2;
      if (k2 > INT_MAX) return INT_MAX;
      k = ceil(max(k1, k2));
    }
    else { // Otherwise, we don't need k2
      k = ceil(k1);
    }

    out[i] = 2 * k + 1; // Convert *pairs* of terms to total *individual* terms
  }
  return out;
}


////////// Navarro and Fuss 2009 ///////////////////////////////////////////////
// [[Rcpp::export]]
IntegerVector ks_Nav(const NumericVector& t, const double& w, const double& eps)
{ // note: w is not used
  int n = t.length();
  IntegerVector out(n);
  double tt;
  float ks, bc;

  for (int i = 0; i < n; i++) {
    tt = t[i];

    if (eps < 1 / (2 * sqrt(2 * M_PI * tt))) { // if error threshold is set low enough
      ks = 2 + sqrt(-2 * tt * log(2 * eps * sqrt(2 * M_PI * tt)));
      bc = sqrt(tt) + 1; // boundary conditions
      if (ks > INT_MAX || bc > INT_MAX) return INT_MAX;
      out[i] = ceil(max(ks, bc)); // ensure boundary conditions are met
    } else {
      out[i] = 2; // else return minimal kappa for that case
    }

  }
  return out;
}
