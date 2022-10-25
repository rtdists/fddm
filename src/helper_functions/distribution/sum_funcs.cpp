// Functions to approximate the infinite sum in the distribution function

#include "declarations.h"



double mills_sum(const double& t, const double& a, const double& v,
                 const double& w, const double& sv, const double& err)
{
  double sum = 0;
  double gamma = v - sv*sv * a * w;
  double lambda = 1 + sv*sv * t;
  double rho = sqrt(t * lambda);
  
  int j = 0;
  double rj = a * j + a * w;
  double m1 = (lambda * rj - gamma * t) / rho;
  double m2 = (lambda * rj + gamma * t) / rho;
  MillsFunc mills_1 = (m1 < 6.5) ? &c_mills : &zeta_mills;
  MillsFunc mills_2 = (m2 < 6.5) ? &c_mills : &zeta_mills;
  double term = SQRT_2PI_INV * exp(-0.5 * rj*rj / t) *
    (mills_1(m1) + mills_2(m2));
  sum += term;
  
  while (term > err) {
    if (j > 1000) {
      warning("pfddm warning: approximation exceeded 1000 terms; the calculation has been stopped and may be inaccurate.");
      break;
    }
    j++;
    rj = a * j + a * (1 - w);
    m1 = (lambda * rj - gamma * t) / rho;
    m2 = (lambda * rj + gamma * t) / rho;
    mills_1 = (m1 < 6.5) ? &c_mills : &zeta_mills;
    mills_2 = (m2 < 6.5) ? &c_mills : &zeta_mills;
    term = SQRT_2PI_INV * exp(-0.5 * rj*rj / t) * (mills_1(m1) + mills_2(m2));
    sum -= term;
    if (term <= err) break;
    j++;
    rj = a * j + a * w;
    m1 = (lambda * rj - gamma * t) / rho;
    m2 = (lambda * rj + gamma * t) / rho;
    mills_1 = (m1 < 6.5) ? &c_mills : &zeta_mills;
    mills_2 = (m2 < 6.5) ? &c_mills : &zeta_mills;
    term = SQRT_2PI_INV * exp(-0.5 * rj*rj / t) * (mills_1(m1) + mills_2(m2));
    sum += term;
  }
  
  return (sum > 0) ? sum : 0; // if result is negative, return 0 instead
}



double ncdf_sum(const double& t, const double& a, const double& v,
                const double& w, const double& sv, const double& err)
{
  double sum = 0;
  double gamma = v - sv*sv * a * w;
  double lambda = 1 + sv*sv * t;
  double rho = sqrt(t * lambda);
  
  int j = 0;
  double rj = a * j + a * w;
  double mult = exp(0.5 * sv*sv * rj*rj);
  if (!isfinite(mult)) {
    warning("pfddm warning: approximation to be multiplied by infinity; the calculation has been stopped and may be inaccurate.");
    return 0;
  }
  double m1 = (gamma * t - lambda * rj) / rho;
  double m2 = (-gamma * t - lambda * rj) / rho;
  double term = mult *  ( exp(-gamma * rj) * R::pnorm(m1, 0, 1, 1, 0) + 
                          exp(gamma * rj) * R::pnorm(m2, 0, 1, 1, 0) );
  sum += term;
  
  while (term > err) {
    if (j > 1000) {
      warning("pfddm warning: approximation exceeded 1000 terms; the calculation has been stopped and may be inaccurate.");
      break;
    }
    j++;
    rj = a * j + a * (1 - w);
    mult = exp(0.5 * sv*sv * rj*rj);
    if (!isfinite(mult)) {
      warning("pfddm warning: approximation to be multiplied by infinity; the calculation has been stopped and may be inaccurate.");
      break;
    }
    m1 = (gamma * t - lambda * rj) / rho;
    m2 = (-gamma * t - lambda * rj) / rho;
    term = mult *  ( exp(-gamma * rj) * R::pnorm(m1, 0, 1, 1, 0) + 
                     exp(gamma * rj) * R::pnorm(m2, 0, 1, 1, 0) );
    sum -= term;
    if (term <= err) break;
    j++;
    rj = a * j + a * w;
    mult = exp(0.5 * sv*sv * rj*rj);
    if (!isfinite(mult)) {
      warning("pfddm warning: approximation to be multiplied by infinity; the calculation has been stopped and may be inaccurate.");
      break;
    }
    m1 = (gamma * t - lambda * rj) / rho;
    m2 = (-gamma * t - lambda * rj) / rho;
    term = mult *  ( exp(-gamma * rj) * R::pnorm(m1, 0, 1, 1, 0) + 
                     exp(gamma * rj) * R::pnorm(m2, 0, 1, 1, 0) );
    sum += term;
  }
  
  return (sum > 0) ? sum : 0; // if result is negative, return 0 instead
}
