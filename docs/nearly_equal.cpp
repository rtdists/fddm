// (Rcpp) C implementation of equality of R vectors

#include <Rcpp.h>
// [[Rcpp::plugins(cpp11)]]
#include <math.h>




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
