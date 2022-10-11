// Helper functions for the distribution function pfddm



void determine_method(const string& method, DisFunc& disf,
                      double& rt0, const bool& log_prob)
{
  char method0 = (!method.empty()) ? method[0] : EMPTYCHAR;
  
  if (log_prob) { // calculate log(probability)
    rt0 = -std::numeric_limits<double>::infinity();
    if (method0 == 'M' || method0 == 'm' || method0 == '1') {
      disf = &Fs_mills_log;
    } else if (method0 == 'N' || method0 == 'n' || method0 == '2') {
      disf = &Fs_ncdf_log;
    } else {
      stop("pfddm error: invalid function parameter 'method': %s.", method);
    }
  } else { // calculate regular probability
    rt0 = 0;
    if (method0 == 'M' || method0 == 'm' || method0 == '1') {
      disf = &Fs_mills;
    } else if (method0 == 'N' || method0 == 'n' || method0 == '2') {
      disf = &Fs_ncdf;
    } else {
      stop("pfddm error: invalid function parameter 'method': %s.", method);
    }
  }
}



double prob_lower(const double& a, const double& v, const double& w,
                  const double& rt0)
{
  double v_threshold = 0.001;
  double prob;
  
  if (rt0 < 0) { // log version, i.e., rt0 = -Inf
    if (-v_threshold < v && v < v_threshold) {
      prob = log(1 - w);
    } else if (v > 0) {
      prob = log( 1 - exp(-2 * v * a * (1 - w)) ) -
             log( exp(2*v*a*w) - exp(-2 * v * a * (1 - w)) );
    } else {
      prob = log( exp(-2 * v * a * (1 - w)) - 1 ) -
             log( exp(-2 * v * a * (1 - w)) - exp(2*v*a*w) );
    }
  } else { // non-log version, i.e., rt0 = 0
    prob = (v < v_threshold && v > -v_threshold) ? 1 - w :
             ( 1 - exp(-2 * v * a * (1 - w)) ) /
             ( exp(2*v*a*w) - exp(-2 * v * a * (1 - w)) );
  }
  
  return prob;
}



void calculate_cdf(const int& Nrt, const int& Nv, const int& Na, const int& Nt0,
                   const int& Nw, const int& Nsv, const int& Nsig,
                   const int& Nerr, const int& Nmax,
                   const NumericVector& rt,
                   const NumericVector& v, const NumericVector& a,
                   const NumericVector& t0, const NumericVector& w,
                   const NumericVector& sv, const NumericVector& sigma,
                   const NumericVector& err, vector<double>& out,
                   const double& rt0, const DisFunc& disf)
{
  double t;
  if (Nsig == 1 && sigma[0] == 1) { // standard diffusion coefficient
    for (int i = 0; i < Nmax; i++) {
      if (isnormal(out[i])) { // not {NaN, NA, Inf, -Inf, rt0 = {0 or -Inf} }
        t = rt[i % Nrt] - t0[i % Nt0]; // response time minus non-decision time
        if (t > 0) { // sort response and calculate density
          if (t > 32) { // appx for +Infinity
            t = 32; // see zedonked/scratch_prob_upper.R for reason it's 32
          }
          if (out[i] == 1) { // response is "lower" so use unchanged parameters
            out[i] = disf(t, v[i % Nv], a[i % Na], w[i % Nw], sv[i % Nsv],
                          err[i % Nerr]);
          } else { // response is "upper" so use alternate parameters
            out[i] = disf(t, -v[i % Nv], a[i % Na], 1 - w[i % Nw],
                          sv[i % Nsv], err[i % Nerr]);
          }
        } else { // {NaN, NA} evaluate to FALSE
          if (isnan(t)) {
            out[i] = t;
          } else {
            out[i] = rt0;
          }
        }
      }
    }
  } else { // non-standard diffusion coefficient
    for (int i = 0; i < Nmax; i++) {
      if (isnormal(out[i])) { // not {NaN, NA, Inf, -Inf, rt0 = {0 or -Inf} }
        t = rt[i % Nrt] - t0[i % Nt0]; // response time minus non-decision time
        if (t > 0) { // sort response and calculate density
          if (t > 32) { // appx for +Infinity
            t = 32; // see zedonked/scratch_prob_upper.R for reason it's 32
          }
          if (out[i] == 1) { // response is "lower" so use unchanged parameters
            out[i] = disf(t, v[i % Nv]/sigma[i % Nsig],
                          a[i % Na]/sigma[i % Nsig], w[i % Nw],
                          sv[i % Nsv]/sigma[i % Nsig], err[i % Nerr]);
          } else { // response is "upper" so use alternate parameters
            out[i] = disf(t, -v[i % Nv]/sigma[i % Nsig],
                          a[i % Na]/sigma[i % Nsig], 1 - w[i % Nw],
                          sv[i % Nsv]/sigma[i % Nsig], err[i % Nerr]);
          }
        } else { // {NaN, NA} evaluate to FALSE
          if (isnan(t)) {
            out[i] = t;
          } else {
            out[i] = rt0;
          }
        }
      }
    }
  }
}
