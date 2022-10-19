// Function to calculate the partial derivative of the DDM PDF

vector<double> partial_pdf(const ParFunc& parf,
                           const NumericVector& rt,
                           const SEXP& response,
                           const NumericVector& v,
                           const NumericVector& a,
                           const NumericVector& t0,
                           const NumericVector& w,
                           const NumericVector& sv,
                           const NumericVector& sigma,
                           const double& sl_thresh,
                           NumericVector err_tol)
{
  // determine lengths of parameter inputs, except response
  int Nrt  = rt.length();
  int Nv   = v.length();
  int Na   = a.length();
  int Nt0  = t0.length();
  int Nw   = w.length();
  int Nsv  = sv.length();
  int Nsig = sigma.length();
  int Nerr = err_tol.length();
  int Nmax = max({Nrt, Nv, Na, Nt0, Nw, Nsv, Nsig, Nerr}); // include Nres later
  int Nres;

  // initialize output, resized in convert_responses() inside parameter_check()
  vector<double> out;

  // check for invalid inputs, invalid inputs get marked in the vector `out`
  double rt0 = 0;
  if (!parameter_check(Nrt, Nres, Nv, Na, Nt0, Nw, Nsv, Nsig, Nerr, Nmax,
                       rt, response, v, a, t0, w, sv, sigma, err_tol,
                       out, rt0)) {
    vector<double> empty_out(0);
    return empty_out;
  }

  // loop through all inputs, the vector `out` gets updated
  double t;
  if (Nsig == 1 && sigma[0] == 1) { // standard diffusion coefficient
    for (int i = 0; i < Nmax; i++) {
      if (isnormal(out[i])) { // not {NaN, NA, Inf, -Inf, rt0 = {0 or -Inf} }
        t = rt[i % Nrt] - t0[i % Nt0]; // response time minus non-decision time
        if (t > 0 && isfinite(t)) {
          if (out[i] == 1) { // response is "lower" so use unchanged parameters
            out[i] = parf(t, v[i % Nv], a[i % Na], w[i % Nw],
                          sv[i % Nsv], err_tol[i % Nerr], sl_thresh);
          } else { // response is "upper" so use alternate parameters
            out[i] = parf(t, -v[i % Nv], a[i % Na], 1 - w[i % Nw],
                          sv[i % Nsv], err_tol[i % Nerr], sl_thresh);
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
              if (t > 0 && isfinite(t)) {
                  if (out[i] == 1) { // response is "lower" so use unchanged parameters
                      out[i] = parf(t, v[i % Nv]/sigma[i % Nsig],
                                    a[i % Na]/sigma[i % Nsig], w[i % Nw],
                                    sv[i % Nsv]/sigma[i % Nsig],
                                    err_tol[i % Nerr], sl_thresh);
                  } else { // response is "upper" so use alternate parameters
                      out[i] = parf(t, -v[i % Nv]/sigma[i % Nsig],
                                    a[i % Na]/sigma[i % Nsig], 1 - w[i % Nw],
                                    sv[i % Nsv]/sigma[i % Nsig],
                                    err_tol[i % Nerr], sl_thresh);
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

  return out;
}
