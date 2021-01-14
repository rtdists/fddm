// Functions to determine the number of terms required for approximating the
// infinite sum in the density function. Note that these functions require the
// scaled time t := t / (a * a), which is input in the function calls in the
// density functions in src/density_funcs.cpp.

#include "funcs.h"




//////////                                                           //////////
///////////////////////////////// Small Time //////////////////////////////////
//////////                                                           //////////


// Gondan et al 2014
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
  return k;
}


// Navarro and Fuss 2009
int ks_Nav(const double& t, const double& w, const double& eps)
{
  // note: w is not used
  if (eps < 1 / (2 * sqrt(2 * M_PI * t))) { // if error threshold is set low enough
    float ks = 2 + sqrt(-2 * t * log(2 * eps * sqrt(2 * M_PI * t)));
    float bc = sqrt(t) + 1; // boundary conditions
    if (ks > INT_MAX || bc > INT_MAX) return INT_MAX;
    return ceil((max(ks, bc) - 1) / 2); // ensure boundary conditions are met
  }
  return 2; // else return minimal kappa for that case
}




//////////                                                           //////////
///////////////////////////////// Large Time //////////////////////////////////
//////////                                                           //////////


// Navarro and Fuss 2009
int kl_Nav(const double& t, const double& w, const double& eps)
{
  // note: w is not used
  float bc = 1 / (M_PI * sqrt(t)); // boundary conditions
  if (bc > INT_MAX) return INT_MAX;
  if (eps < 1 / (M_PI * t)) { // error threshold is low enough
    float kl = sqrt(-2 * log(M_PI * t * eps) / (M_PI*M_PI * t));
    if (kl > INT_MAX) return INT_MAX;
    return ceil(max(kl, bc)); // ensure boundary conditions are met
  }
  return ceil(bc); // else set to boundary condition
}
