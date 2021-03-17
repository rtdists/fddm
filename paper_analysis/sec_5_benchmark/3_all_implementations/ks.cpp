// C++ functions for calculating the required number of terms in the truncated
// small-time sum. These functions are based on the ones in `src/num_funcs.cpp`,
// but they accept a vectors as inputs and loop over them. This is done to add
// more computation time in C++ so that the `microbenchmark` results contain
// proportionally less overhead from R.

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
IntegerVector ks_Gon(const NumericVector& t, const NumericVector& w,
                     const NumericVector& eps)
{
  int nt = t.length();
  int nw = w.length();
  int ne = eps.length();
  int n = nt * nw * ne;
  IntegerVector out(n);
  float u_eps, arg, k1;
  int k;

  for (int ti = 0; ti < nt; ti++) {
    for (int wi = 0; wi < nw; wi++) {
      for (int ei = 0; ei < ne; ei++) {
        u_eps = min(-1.0, log(2 * M_PI * t[ti]*t[ti] * eps[ei]*eps[ei])); // Safe bound for sqrt
        arg = -t[ti] * (u_eps - sqrt(-2 * u_eps - 2)); // sqrt(x) with x > 0
        k1 = (sqrt(2 * t[ti]) - w[wi])/2;
        if (k1 > INT_MAX) return INT_MAX;
        if (arg > 0) { // If arg > 0, set k2 and calculate k
          float k2 = (sqrt(arg) - w[wi]) / 2;
          if (k2 > INT_MAX) return INT_MAX;
          k = ceil(max(k1, k2));
        }
        else { // Otherwise, we don't need k2
          k = ceil(k1);
        }

        // Convert *pairs* of terms to total *individual* terms
        out[ti*nw*ne + wi*ne + ei] = 2 * k + 1;
      }
    }
  }
  return out;
}


////////// Navarro and Fuss 2009 ///////////////////////////////////////////////
// [[Rcpp::export]]
IntegerVector ks_Nav(const NumericVector& t, const NumericVector& w,
                     const NumericVector& eps)
{ // note: w is not used
  int nt = t.length();
  int nw = w.length();
  int ne = eps.length();
  int n = nt * nw * ne;
  IntegerVector out(n);
  float ks, bc;

  for (int ti = 0; ti < nt; ti++) {
    for (int wi = 0; wi < nw; wi++) {
      for (int ei = 0; ei < ne; ei++) {
        if (eps[ei] < 1 / (2 * sqrt(2 * M_PI * t[ti]))) { // if error threshold is set low enough
          ks = 2 + sqrt(-2 * t[ti] * log(2 * eps[ei] * sqrt(2 * M_PI * t[ti])));
          bc = sqrt(t[ti]) + 1; // boundary conditions
          if (ks > INT_MAX || bc > INT_MAX) return INT_MAX;
          out[ti*nw*ne + wi*ne + ei] = ceil(max(ks, bc)); // ensure boundary conditions are met
        } else {
          out[ti*nw*ne + wi*ne + ei] = 2; // else return minimal kappa for that case
        }
      }
    }
  }
  return out;
}
