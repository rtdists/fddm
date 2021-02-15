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


////////// fddm 2021 ///////////////////////////////////////////////////////////
// [[Rcpp::export]]
IntegerVector ks_SWSE_14(const NumericVector& t, const NumericVector& w,
                         const NumericVector& eps)
{ // note: t is assumed to be t := t / (a*a)
  int nt = t.length();
  int nw = w.length();
  int ne = eps.length();
  IntegerVector out(nt * nw * ne);

  double gamma, sum, term, tt, ww, eeps, rj;
  int minterms, j, oddj;

  for (int ti = 0; ti < nt; ti++) {
    tt = t[ti];
    for (int wi = 0; wi < nw; wi++) {
      ww = w[wi];
      for (int ei = 0; ei < ne; ei++) {
        eeps = eps[ei];

        minterms = sqrt(tt)/2 - ww/2; // min number of terms, truncates toward 0
        gamma = -1 / (2 * tt);
        sum = ww * exp(gamma * ww*ww); // initialize with j=0 term
        j = 0;
        while (j < minterms) {
          j++;
          rj = 2*j - ww;
          sum -= rj * exp(gamma * rj*rj);
          rj = 2*j + ww;
          sum += rj * exp(gamma * rj*rj);
        }
        j++;
        rj = 2*j - ww;
        term = rj * exp(gamma * rj*rj);
        sum -= term;
        oddj = 1;
        while (term > eeps) {
          rj = 2*j + ww;
          term = rj * exp(gamma * rj*rj);
          sum += term;
          if (term <= eeps) {
            oddj = 0;
            break;
          }
          j++;
          rj = 2*j - ww;
          term = rj * exp(gamma * rj*rj);
          sum -= term;
        }
        out[ti*nw*ne + wi*ne + ei] = 2*j - oddj;
      }
    }
  }
  return out;
}

// [[Rcpp::export]]
IntegerVector ks_SWSE_17(const NumericVector& t, const NumericVector& w,
                         const NumericVector& eps)
{ // note: t is assumed to be t := t / (a*a)
  int nt = t.length();
  int nw = w.length();
  int ne = eps.length();
  IntegerVector out(nt * nw * ne);

  double gamma, sum, term, rj, tt, ww, eeps;
  int minterms, j;

  for (int ti = 0; ti < nt; ti++) {
    tt = t[ti];
    for (int wi = 0; wi < nw; wi++) {
      ww = w[wi];
      for (int ei = 0; ei < ne; ei++) {
        eeps = eps[ei];

        minterms = sqrt(tt) - ww; // min number of terms, truncates toward 0
        gamma = -1 / (2 * tt);
        sum = ww * exp(gamma * ww*ww); // initialize with j=0 term
        j = 0;
        if (minterms % 2) { // minterms is odd (and at least 1)
          j++;
          rj = j + 1 - ww;
          term = rj * exp(gamma * rj*rj);
          sum -= term;
          while (j < minterms) {
            j++;
            rj = j + ww;
            sum += rj * exp(gamma * rj*rj);
            j++;
            rj = j + 1 - ww;
            term = rj * exp(gamma * rj*rj);
            sum -= term;
          }
          j++;
          rj = j + ww; // j is now even
          term = rj * exp(gamma * rj*rj);
          sum += term;
          while (term > eeps) {
            j++;
            rj = j + 1 - ww;
            term = rj * exp(gamma * rj*rj);
            sum -= term;
            if (term <= eeps) break;
            j++;
            rj = j + ww;
            term = rj * exp(gamma * rj*rj);
            sum += term;
          }
        } else { // minterms is even (and at least 0)
          while (j < minterms) { // j is currently 0
            j++;
            rj = j + 1 - ww;
            sum -= rj * exp(gamma * rj*rj);
            j++;
            rj = j + ww;
            term = rj * exp(gamma * rj*rj);
            sum += term;
          }
          j++;
          rj = j + 1 - ww; // j is now odd
          term = rj * exp(gamma * rj*rj);
          sum -= term;
          while (term > eeps) {
            j++;
            rj = j + ww;
            term = rj * exp(gamma * rj*rj);
            sum += term;
            if (term <= eeps) break;
            j++;
            rj = j + 1 - ww;
            term = rj * exp(gamma * rj*rj);
            sum -= term;
          }
        }
        out[ti*nw*ne + wi*ne + ei] = j;
      }
    }
  }
  return out;
}
