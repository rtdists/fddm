// Functions to determine the number of terms required for approximating the
// infinite sum in the density function. Note that these functions require the
// scaled time t := t / (a * a), which is input in the function calls in the
// density functions in src/density_funcs.cpp.



///////////////////////////////// Large Time ///////////////////////////////////
// number of terms required for the large-time sum, not differentiated
// (Navarro and Fuss 2009)
int kl_Nav(const double& taa, const double& err)
{ // note: taa = t / (a*a)
  float bc = O_PI / sqrt(taa); // boundary conditions
  if (bc > INT_MAX) return INT_MAX;
  if (err * PI_CONST * taa < 1) { // error threshold is low enough
    float kl = sqrt(-2 * log(PI_CONST * taa * err) / (PI2 * taa));
    if (kl > INT_MAX) return INT_MAX;
    return ceil(max(kl, bc)); // ensure boundary conditions are met
  }
  return ceil(bc); // else set to boundary condition
}

// number of terms required for the large-time sum, differentiated
// (Hartmann and Klauer 2021)
int kl_Har(const double& taa, const double& t, const double& err)
{ // note: taa = t / (a*a)
  float k1 = SQRT_3 * O_PI / sqrt(taa);
  double ueps = log(0.6 * taa * t * PI_CONST * err);
  float k2 = SQRT_2_PI * sqrt((-ueps + SQRT_2 * sqrt(-ueps - 1)) / taa);
  if (k1 > INT_MAX || k2 > INT_MAX) {
    return INT_MAX;
  } else {
    return ceil(max(k1, k2));
  }
}

// number of terms required for the large-time sum for w, differentiated
// (Hartmann and Klauer 2021)
int kl_Har_w(const double& taa, const double& t, const double& err)
{ // note: taa = t / (a*a)
  float k1 = SQRT_2_1_PI / taa;
  double ueps = log(FOUR_NINTHS * taa*taa*taa * PI_CONST*PI_CONST * err*err);
  float k2 = O_PI * sqrt((-ueps + SQRT_2 * sqrt(-ueps - 1)) / sqrt(taa));
  if (k1 > INT_MAX || k2 > INT_MAX) {
    return INT_MAX;
  } else {
    return ceil(max(k1, k2));
  }
}



///////////////////////////////// Small Time ///////////////////////////////////
// number of terms required for the small-time sum for w, differentiated
// (Hartmann and Klauer 2021)
int ks_Har_w(const double& taa, const double& w, const double& err)
{
  float k1 = 0.5 * w + 0.5 * SQRT_3 * sqrt(taa);
  double ueps = -LOG_2 - LOG_PI - 2 * log(err) - 2 * log(taa);
  float k2 = 0.5 * w + 0.5 * sqrt(taa * (ueps + SQRT_2 * sqrt(ueps - 1)));
  if (k1 > INT_MAX || k2 > INT_MAX) {
    return INT_MAX;
  } else {
    return ceil(max(k1, k2));
  }
}
