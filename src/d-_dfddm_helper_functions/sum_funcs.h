// Functions to approximate the infinite sum in the derivative of the PDF



///////////////////////////////// Small Time ///////////////////////////////////
// sum for small time (not differentiated)
// adapted from (Blurton et al 2017)
double sum_small(const double& taa, const double& w, const double& err)
{ // note: taa = t / (a*a)
  int minterms = sqrt(taa) - w; // min number of terms, truncates toward 0
  double gamma = -0.5 / taa;
  double sum = w * exp(gamma * w*w); // initialize with j=0 term
  double term, rj;
  int j = 0;
  if (minterms % 2) { // minterms is odd (and at least 1)
    j++; // j is now 1
    rj = j + 1 - w;
    sum -= rj * exp(gamma * rj*rj);
    while (j < minterms) {
      j++;
      rj = j + w;
      sum += rj * exp(gamma * rj*rj);
      j++;
      rj = j + 1 - w;
      sum -= rj * exp(gamma * rj*rj);
    }
    j++; // need at least the next term to check if small enough
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
      sum += rj * exp(gamma * rj*rj);
    }
    j++; // need at least the next term to check if small enough
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
  return sum;
}

// sum for small time (differentiated)
// based on (Blurton et al 2017)
double sum_small_d(const double& taa, const double& w, const double& err)
{ // note: taa = t / (a*a)
  int minterms = SQRT_3 * sqrt(taa) - w; // min number of terms, truncates toward 0
  double gamma = -0.5 / taa;
  double sum = w*w*w * exp(gamma * w*w); // initialize with j=0 term
  double term, rj;
  int j = 0;
  if (minterms % 2) { // minterms is odd (and at least 1)
    j++; // j is now 1
    rj = j + 1 - w;
    sum -= rj*rj*rj * exp(gamma * rj*rj);
    while (j < minterms) {
      j++;
      rj = j + w;
      sum += rj*rj*rj * exp(gamma * rj*rj);
      j++;
      rj = j + 1 - w;
      sum -= rj*rj*rj * exp(gamma * rj*rj);;
    }
    j++; // need at least the next term to check if small enough
    rj = j + w; // j is now even
    term = rj*rj*rj * exp(gamma * rj*rj);
    sum += term;
    while (term > err) {
      j++;
      rj = j + 1 - w;
      term = rj*rj*rj * exp(gamma * rj*rj);
      sum -= term;
      if (term <= err) break;
      j++;
      rj = j + w;
      term = rj*rj*rj * exp(gamma * rj*rj);
      sum += term;
    }
  } else { // minterms is even (and at least 0)
    while (j < minterms) { // j is currently 0
      j++;
      rj = j + 1 - w;
      sum -= rj*rj*rj * exp(gamma * rj*rj);
      j++;
      rj = j + w;
      sum += rj*rj*rj * exp(gamma * rj*rj);
    }
    j++; // need at least the next term to check if small enough
    rj = j + 1 - w; // j is now odd
    term = rj*rj*rj * exp(gamma * rj*rj);
    sum -= term;
    while (term > err) {
      j++;
      rj = j + w;
      term = rj*rj*rj * exp(gamma * rj*rj);
      sum += term;
      if (term <= err) break;
      j++;
      rj = j + 1 - w;
      term = rj*rj*rj * exp(gamma * rj*rj);
      sum -= term;
    }
  }
  return sum;
}

// sum for small time for dw (differentiated)
// adapted from (Hartmann and Klauer 2021)
double sum_small_d_w(const double& taa, const double& w, const int& ks)
{
  double gamma = -1 / taa;
  double rj = gamma * w*w;
  double sum = (1 + rj) * exp(0.5 * rj); // initialize with j = 0 term
  for (int j = 1; j <= ks; j++) {
    rj = gamma * (w + 2 * j)*(w + 2 * j);
    sum += (1 + rj) * exp(0.5 * rj);
    rj = gamma * (w - 2 * j)*(w - 2 * j);
    sum += (1 + rj) * exp(0.5 * rj);
  }
  return sum;
}


// sum for log-differentiated
// double sum_small_log_d(const double& taa, const double& w, const double& err)
// { // note: taa = t / (a*a)
  
// }


///////////////////////////////// Large Time ///////////////////////////////////
// sum for large time (not differentiated)
// (Navarro and Fuss 2009)
double sum_large(const double& taa, const double& w, const int& kl)
{ // note: taa = t / (a*a)
  double gamma = -0.5 * PI_CONST*PI_CONST * taa;
  double sum = 0.0;
  for (int j = 1; j <= kl; j++) {
    sum += j * sin(j * w * PI_CONST) * exp(gamma * j*j);
  }
  return sum;
}

// sum for large time (differentiated)
// adapted from (Hartmann and Klauer 2021)
double sum_large_d(const double& taa, const double& w, const int& kl)
{ // note: taa = t / (a*a)
  double gamma = -0.5 * PI2 * taa;
  double sum = 0.0;
  for (int j = 1; j <= kl; j++) {
    sum += j*j*j * sin(j * w * PI_CONST) * exp(gamma * j*j);
  }
  return sum;
}

// sum for large time for dw (differentiated)
// adapted from (Hartmann and Klauer 2021)
double sum_large_d_w(const double& taa, const double& w, const int& kl)
{ // note: taa = t / (a*a)
  double gamma = -0.5 * PI2 * taa;
  double sum = 0.0;
  for (int j = 1; j <= kl; j++) {
    sum += j*j * cos(j * w * PI_CONST) * exp(gamma * j*j);
  }
  return sum;
}
