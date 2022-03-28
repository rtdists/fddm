// Functions to calculate the number of terms required in the infinite sums
// included from src/fitting_helper_functions/likelihood_funcs.h
// included from src/fitting_helper_functions/gradient_funcs.h



//--------------------- Large-Time -------------------------------------------//
// for Large-Time PDF
int kl_pdf(const double& taa, const double& err)
{ // note: taa = t / (a*a)
  double bc = 1 / (PI_CONST * sqrt(taa)); // boundary conditions
  if (bc > INT_MAX) return INT_MAX;
  if (err * PI_CONST * taa < 1) { // error threshold is low enough
    double kl = sqrt(-2 * log(PI_CONST * taa * err)
                     / (PI_CONST*PI_CONST * taa));
    if (kl > INT_MAX) {
      return INT_MAX;
    }
    return ceil(max(kl, bc)); // ensure boundary conditions are met
  }
  return ceil(bc); // else set to boundary condition
}

// for Large-Time derivative wrt a, t0
int kl_dat(const double& taa, const double& t, const double& err)
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

// for Large-Time derivative wrt w
int kl_dw(const double& taa, const double& t, const double& err)
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



//--------------------- Small-Time -------------------------------------------//
// for Small-Time derivative wrt w
int ks_dw(const double& taa, const double& w, const double& err)
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
