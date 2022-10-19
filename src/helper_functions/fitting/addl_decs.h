#include "../parameter_checks.h"

using std::isnormal;
using Rcpp::NumericVector;

//--------- Templates for Partial Derivatives --------------------------------//
typedef double (*ParFunc)(const double&, const double&, const double&,
                          const double&, const double&, const double&,
                          const double&);
