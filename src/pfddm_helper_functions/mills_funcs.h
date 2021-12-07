// Functions to calculate (or approximate) Mills ratio



double r_mills(const double& x)
{
  return R::pnorm(x, 0, 1, 0, 0) / R::dnorm(x, 0, 1, 0);
}

double c_mills(const double& x)
{
  return SQRT_2PI * 0.5 * (1 + erf(SQRT_2_INV_NEG * x)) * exp(0.5 * x*x);;
}

double zeta_mills(const double& x)
{
  double x2 = x*x;
  return ( 1 -
           1 / (x2 + 2) +
           1 / ( (x2 + 2) * (x2 + 4) ) -
           5 / ( (x2 + 2) * (x2 + 4) * (x2 + 6) ) +
           9 / ( (x2 + 2) * (x2 + 4) * (x2 + 6) * (x2 + 8) ) -
           129 / ( (x2 + 2) * (x2 + 4) * (x2 + 6) * (x2 + 8) * (x2 + 10) )
  ) / x;
}