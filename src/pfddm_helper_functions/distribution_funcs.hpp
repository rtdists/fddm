// Functions to evaluate the DDM CDF for specific criteria



// 1. with sv Mills ratio
double Fs_mills(const double& t, const double& a, const double& v,
                const double& w, const double& sv, const double& err)
{
  double mult = exp( (sv*sv * a*a * w*w - 2 * v * a * w - v*v * t) /
                     (2 + 2 * sv*sv * t) );
  double out = mult * mills_sum(t, a, v, w, sv, err / mult);
  
  return (out < 1) ? out : 1; // occasionally cdf is slightly larger than 1
}

double Fs_mills_log(const double& t, const double& a, const double& v,
                    const double& w, const double& sv, const double& err)
{
  double mult = (sv*sv * a*a * w*w - 2 * v * a * w - v*v * t) /
    (2 + 2 * sv*sv * t);
  
  double temp = mills_sum(t, a, v, w, sv, err / exp(mult));
  if (temp > 0) {
    return mult + log(temp);
  } else{ // protect against -Inf
    return log(err) - LOG_100;
  }
}



// 2. with sv normal CDF
double Fs_ncdf(const double& t, const double& a, const double& v,
               const double& w, const double& sv, const double& err)
{
  double mult = exp(0.5 * sv*sv * a*a * w*w - v * a * w);
  double out = mult * ncdf_sum(t, a, v, w, sv, err / mult);
  
  return (out < 1) ? out : 1; // occasionally cdf is slightly larger than 1
}

double Fs_ncdf_log(const double& t, const double& a, const double& v,
                   const double& w, const double& sv, const double& err)
{
  double mult = 0.5 * sv*sv * a*a * w*w - v * a * w;
  
  double temp = ncdf_sum(t, a, v, w, sv, err / exp(mult));
  if (temp > 0) {
    return mult + log(temp);
  } else{ // protect against -Inf
    return log(err) - LOG_100;
  }
}