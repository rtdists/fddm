// (Rcpp) C implementation of equality of R vectors

#include <iostream>
#include <Rcpp.h>
// [[Rcpp::plugins(cpp11)]]
#include <math.h>




//
// Find the maximum element-wise difference between two R vectors
//
// [[Rcpp::export]]
double max_error(Rcpp::NumericVector a, Rcpp::NumericVector b, int n=0) {
  if (n == 0) {
    n = a.size(); // Get number of values in a
    if (b.size() != n) { // check length of a against length of b
      std::cout << "Vectors are not the same length" << std::endl;
      return -1.0;
    }
  }

  double max_err = 0.0;
  for (int i = 0; i < n; i++) {
    max_err = std::max(max_err, fabs(a[i]-b[i])); // Get largest difference
  }

  return max_err;
}
// /*** R
// a = seq(from=0, to=1, by=0.05)
// b = seq(from=0.02, to=1.02, by=0.05)
// n = length(a)
// max_error(a=a, b=b, n=n)
// */


//
// Check if the values of two R vectors are approximately equal
//
// [[Rcpp::export]]
bool almost_equal(Rcpp::NumericVector a, Rcpp::NumericVector b,
                  double thresh = 0.0) {
  if (thresh == 0.0) { // set default value (doesn't work in definition)
    thresh = sqrt(DBL_EPSILON);
  }

  int n = a.size(); // Get number of values in a
  if (b.size() != n) { // check length of a against length of b
    std::cout << "Vectors are not the same length" << std::endl;
    return false;
  }

  double max_err = 0.0;
  for (int i = 0; i < n; i++) {
    max_err = std::max(max_err, fabs(a[i]-b[i])); // Get largest difference
  }

  if (max_err <= thresh) {
    return true;
  }
  else {
    std::cout << "max error is: " << max_err << std::endl;
    return false;
  }
}
// /*** R
// a = seq(from=0, to=1, by=0.05)
// b = seq(from=0.02, to=1.02, by=0.05)
// thresh = 0.03
// almost_equal(a=a, b=b, thresh=thresh)
// almost_equal(a=a, b=b)
// */


//
// Check if all the values in an R vector are approximately equal;
//   -treat first element as the "truth"
//
// [[Rcpp::export]]
Rcpp::NumericVector nearly_equal(Rcpp::NumericVector a, double thresh = 0.0) {
  if (thresh == 0.0) { // set default value (doesn't work in definition)
    thresh = sqrt(DBL_EPSILON);
  }
  thresh *= 2; // allow for approach from either side of limit (of "truth")

  int n = a.size(); // Get number of elements in the input vector
  Rcpp::NumericVector out(n); // Initialize output vector
  out[0] = 0; // set value of the first element as the "truth"

  if (n > 1) { // loop through all elements of vector
    for (int i = 1; i < n; i++) {
      double dif = fabs(a[0]-a[i]);
      if (dif <= thresh) {
        out[i] = 0;
      } else {
        out[i] = dif;
      }
    }
  }

  return out;
}
