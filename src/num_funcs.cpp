// Functions to determine the number of terms required for approximating the
// infinite sum in the density function. Note that these functions require the
// scaled time t := t / (a * a), which is input in the function calls in the
// density functions in src/density_funcs.cpp.

#include "funcs.h"




//////////                                                           //////////
///////////////////////////////// Small Time //////////////////////////////////
//////////                                                           //////////


// Kesselmeier and co 2014
int ks_Kes(const double& t, const double& w, const double& eps)
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
  return k;
}


// Navarro2009
int ks_Nav(const double& t, const double& w, const double& eps)
{
  // note: w is not used
  if (eps < 1 / (2 * sqrt(2 * M_PI * t))) { // if error threshold is set low enough
    double ks = 2 + sqrt(-2 * t * log(2 * eps * sqrt(2 * M_PI * t)));
    return ceil(max(ks, sqrt(t)+1)); // ensure boundary conditions are met
  }
  return 2; // else return minimal kappa for that case
}




//////////                                                           //////////
///////////////////////////////// Large Time //////////////////////////////////
//////////                                                           //////////


// Navarro2009
int kl_Nav(const double& t, const double& w, const double& eps)
{
  // note: w is not used
  if (eps < 1 / (M_PI * t)) { // error threshold is low enough
    double kl = sqrt(-2 * log(M_PI * t * eps) / (M_PI*M_PI *t));
    return ceil(max(kl, 1 / (M_PI * sqrt(t))));
  }
  return ceil(1 / (M_PI * sqrt(t))); // else set to boundary condition
}
