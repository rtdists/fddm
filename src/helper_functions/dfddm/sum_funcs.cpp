// Functions to approximate the infinite sum in the density function

#include "declarations.h"


//////////                                                           //////////
///////////////////////////////// Small Time //////////////////////////////////
//////////                                                           //////////


// term < err Blurton et al 2017 style truncated sum, with minimum terms
double small_sum_eps_17(const double& t, const double& a, const double& w,
                        const int& ks, const double& err)
{ // note: ks is not used
  int minterms = sqrt(t)/a - w; // min number of terms, truncates toward 0
  double gamma = -a*a / (2 * t);
  double sum = w * exp(gamma * w*w); // initialize with j=0 term
  double term, rj;
  int j = 0;
  if (minterms % 2) { // minterms is odd (and at least 1)
    j++;
    rj = j + 1 - w;
    term = rj * exp(gamma * rj*rj);
    sum -= term;
    while (j < minterms) {
      j++;
      rj = j + w;
      sum += rj * exp(gamma * rj*rj);
      j++;
      rj = j + 1 - w;
      term = rj * exp(gamma * rj*rj);
      sum -= term;
    }
    j++;
    rj = j + w; // j is now even
    term = rj * exp(gamma * rj*rj);
    sum += term;
    while (term > err) {
      j++;
      rj = j + 1 - w;
      term = rj * exp(gamma * rj*rj);
      sum -= term;
      if (term <= err) break;
      j++;
      rj = j + w;
      term = rj * exp(gamma * rj*rj);
      sum += term;
    }
  } else { // minterms is even (and at least 0)
    while (j < minterms) { // j is currently 0
      j++;
      rj = j + 1 - w;
      sum -= rj * exp(gamma * rj*rj);
      j++;
      rj = j + w;
      term = rj * exp(gamma * rj*rj);
      sum += term;
    }
    j++;
    rj = j + 1 - w; // j is now odd
    term = rj * exp(gamma * rj*rj);
    sum -= term;
    while (term > err) {
      j++;
      rj = j + w;
      term = rj * exp(gamma * rj*rj);
      sum += term;
      if (term <= err) break;
      j++;
      rj = j + 1 - w;
      term = rj * exp(gamma * rj*rj);
      sum -= term;
    }
  }
  return (sum > 0) ? sum : 0; // if result is negative, return 0 instead
}


// term < err Gondan et al 2014 style truncated sum, with minimum terms
double small_sum_eps_14(const double& t, const double& a, const double& w,
                        const int& ks, const double& err)
{ // note: ks is not used
  int minterms = sqrt(t)/a/2 - w/2; // min number of terms, truncates toward 0
  double gamma = -a*a / (2 * t);
  double sum = w * exp(gamma * w*w); // initialize with j=0 term
  double term, rj;
  int j = 0;
  while (j < minterms) {
    j++;
    rj = 2*j - w;
    sum -= rj * exp(gamma * rj*rj);
    rj = 2*j + w;
    sum += rj * exp(gamma * rj*rj);
  }
  j++;
  rj = 2*j - w;
  term = rj * exp(gamma * rj*rj);
  sum -= term;
  while (term > err) {
    rj = 2*j + w;
    term = rj * exp(gamma * rj*rj);
    sum += term;
    if (term <= err) {
      break;
    }
    j++;
    rj = 2*j - w;
    term = rj * exp(gamma * rj*rj);
    sum -= term;
  }
  return (sum > 0) ? sum : 0; // if result is negative, return 0 instead
}


// Blurton et al 2017 style truncated sum
double small_sum_2017(const double& t, const double& a, const double& w,
                      const int& ks, const double& err)
{ // note: err is not used
  double gamma = -a*a / (2 * t);
  double sum = w * exp(gamma * w*w); // initialize with j = 0 term
  double rj;
  int j = (ks > 1) ? ks - 1: 0; // check ks
  if (j % 2) { // j is odd (i.e., ks - 1 is odd)
    rj = j + 1 - w;
    sum -= rj * exp(gamma * rj*rj);
    j--;
  }
  while (j > 0) { // iterate from j = ks-1 to j = 1
    rj = j + w; // j is even
    sum += rj * exp(gamma * rj*rj);
    j--;
    rj = j + 1 - w; // j is odd
    sum -= rj * exp(gamma * rj*rj);
    j--;
  }
  return (sum > 0) ? sum : 0; // if result is negative, return 0 instead
}


// Gondan et al 2014 style truncated sum
double small_sum_2014(const double& t, const double& a, const double& w,
                      const int& ks, const double& err)
{ // note: err is not used
  double gamma = -a*a / (2 * t);
  double sum = w * exp(gamma * w*w); // initialize with j = 0 term
  for (int j = ks/2; j > 0; j--) { // iterate over ks/2, which casts to int and
                                   // truncates toward 0 (same as floor(ks/2.0))
    sum += (2 * j + w) * exp(gamma * (2 * j + w) * (2 * j + w))
    - (2 * j - w) * exp(gamma * (2 * j - w) * (2 * j - w));
  }
  return (sum > 0) ? sum : 0; // if result is negative, return 0 instead
}



//////////                                                           //////////
///////////////////////////////// Large Time //////////////////////////////////
//////////                                                           //////////


// Navarro and Fuss 2009 style truncated sum
double large_sum_Nav(const double& t, const double& a, const double& w,
                     const int& kl, const double& err)
{ // note: err is not used
  double gamma = -PI_CONST*PI_CONST * t / (2 * a*a);
  double sum = 0.0;
  for (int j = 1; j <= kl; j++) {
    sum += sin(j * w * PI_CONST) * exp(gamma * j*j) * j;
  }
  return (sum > 0) ? sum : 0; // if result is negative, return 0 instead
}
