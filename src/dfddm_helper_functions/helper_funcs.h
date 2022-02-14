// Helper functions for the density function dfddm



void determine_method(const string& n_terms_small,
                      const string& summation_small,
                      const string& switch_mech, double& switch_thresh,
                      NumFunc& numf, SumFunc& sumf, DenFunc& denf,
                      double& rt0, const bool& log_prob)
{
  if (log_prob) { // calculate log(density)
    rt0 = -std::numeric_limits<double>::infinity();
    if (switch_mech == "eff_rt") {
      if (switch_thresh < 0) switch_thresh = 0.8; // default for this mechanism
      denf = &ft_log;
      if (summation_small == "2017") {
        sumf = &small_sum_eps_17;
      } else if (summation_small == "2014") {
        sumf = &small_sum_eps_14;
      } else {
        stop("dfddm error: invalid function parameter 'summation_small': %s.",
             summation_small);
      }
    } else if (switch_mech == "terms_large") {
      if (switch_thresh < 0) switch_thresh = 1; // default for this mechanism
      denf = &fc_log;
      if (summation_small == "2017") {
        sumf = &small_sum_eps_17;
      } else if (summation_small == "2014") {
        sumf = &small_sum_eps_14;
      } else {
        stop("dfddm error: invalid function parameter 'summation_small': %s.",
             summation_small);
      }
    } else if (switch_mech == "terms") {
      denf = &fb_log;
      if (n_terms_small == "Gondan") {
        numf = &ks_Gon;
      } else if (n_terms_small == "Navarro") {
        numf = &ks_Nav;
      } else {
        stop("dfddm error: invalid function parameter 'n_terms_small': %s.",
             n_terms_small);
      }
      if (summation_small == "2017") {
        sumf = &small_sum_2017;
      } else if (summation_small == "2014") {
        sumf = &small_sum_2014;
      } else {
        stop("dfddm error: invalid function parameter 'summation_small': %s.",
             summation_small);
      }
    } else if (switch_mech == "small") {
      if (n_terms_small == "SWSE") {
        denf = &ff_log;
        if (summation_small == "2017") {
          sumf = &small_sum_eps_17;
        } else if (summation_small == "2014") {
          sumf = &small_sum_eps_14;
        } else {
          stop("dfddm error: invalid function parameter 'summation_small': %s.",
              summation_small);
        }
      } else {
        denf = &fs_log;
        if (n_terms_small == "Gondan") {
          numf = &ks_Gon;
        } else if (n_terms_small == "Navarro") {
          numf = &ks_Nav;
        } else {
          stop("dfddm error: invalid function parameter 'n_terms_small': %s.",
             n_terms_small);
        }
        if (summation_small == "2017") {
          sumf = &small_sum_2017;
        } else if (summation_small == "2014") {
          sumf = &small_sum_2014;
        } else {
          stop("dfddm error: invalid function parameter 'summation_small': %s.",
              summation_small);
        }
      }
    } else if (switch_mech == "large") {
      denf = &fl_log;
    } else {
      stop("dfddm error: invalid function parameter 'switch_mech': %s.",
           switch_mech);
    }
  } else { // calculate regular (non-log) density
    rt0 = 0;
    if (switch_mech == "eff_rt") {
      if (switch_thresh < 0) switch_thresh = 0.8; // default for this mechanism
      denf = &ft;
      if (summation_small == "2017") {
        sumf = &small_sum_eps_17;
      } else if (summation_small == "2014") {
        sumf = &small_sum_eps_14;
      } else {
        stop("dfddm error: invalid function parameter 'summation_small': %s.",
             summation_small);
      }
    } else if (switch_mech == "terms_large") {
      if (switch_thresh < 0) switch_thresh = 1; // default for this mechanism
      denf = &fc;
      if (summation_small == "2017") {
        sumf = &small_sum_eps_17;
      } else if (summation_small == "2014") {
        sumf = &small_sum_eps_14;
      } else {
        stop("dfddm error: invalid function parameter 'summation_small': %s.",
             summation_small);
      }
    } else if (switch_mech == "terms") {
      denf = &fb;
      if (n_terms_small == "Gondan") {
        numf = &ks_Gon;
      } else if (n_terms_small == "Navarro") {
        numf = &ks_Nav;
      } else {
        stop("dfddm error: invalid function parameter 'n_terms_small': %s.",
             n_terms_small);
      }
      if (summation_small == "2017") {
        sumf = &small_sum_2017;
      } else if (summation_small == "2014") {
        sumf = &small_sum_2014;
      } else {
        stop("dfddm error: invalid function parameter 'summation_small': %s.",
             summation_small);
      }
    } else if (switch_mech == "small") {
      if (n_terms_small == "SWSE") {
        denf = &ff;
        if (summation_small == "2017") {
          sumf = &small_sum_eps_17;
        } else if (summation_small == "2014") {
          sumf = &small_sum_eps_14;
        } else {
          stop("dfddm error: invalid function parameter 'summation_small': %s.",
              summation_small);
        }
      } else {
        denf = &fs;
        if (n_terms_small == "Gondan") {
          numf = &ks_Gon;
        } else if (n_terms_small == "Navarro") {
          numf = &ks_Nav;
        } else {
          stop("dfddm error: invalid function parameter 'n_terms_small': %s.",
             n_terms_small);
        }
        if (summation_small == "2017") {
          sumf = &small_sum_2017;
        } else if (summation_small == "2014") {
          sumf = &small_sum_2014;
        } else {
          stop("dfddm error: invalid function parameter 'summation_small': %s.",
              summation_small);
        }
      }
    } else if (switch_mech == "large") {
      denf = &fl;
    } else {
      stop("dfddm error: invalid function parameter 'switch_mech': %s.",
           switch_mech);
    }
  }
}



void calculate_pdf(const int& Nrt, const int& Na, const int& Nv, const int& Nt0,
                   const int& Nw, const int& Nsv, const int& Nsig,
                   const int& Nerr, const int& Nmax,
                   const NumericVector& rt,
                   const NumericVector& a, const NumericVector& v,
                   const NumericVector& t0, const NumericVector& w,
                   const NumericVector& sv, const NumericVector& sigma,
                   const NumericVector& err, vector<double>& out,
                   const double& switch_thresh,
                   const NumFunc& numf, const SumFunc& sumf,
                   const DenFunc& denf, const double& rt0)
{
  double t;
  if (Nsig == 1 && sigma[0] == 1) { // standard diffusion coefficient
    for (int i = 0; i < Nmax; i++) {
      if (isnormal(out[i])) { // not {NaN, NA, Inf, -Inf, rt0 = {0 or -Inf} }
        t = rt[i % Nrt] - t0[i % Nt0]; // response time minus non-decision time
        if (t > 0 && isfinite(t)) { // sort response and calculate density
          if (out[i] == 1) { // response is "lower" so use unchanged parameters
              out[i] = denf(t, a[i % Na], v[i % Nv], w[i % Nw], sv[i % Nsv],
                          err[i % Nerr], switch_thresh, numf, sumf);
          } else { // response is "upper" so use alternate parameters
            out[i] = denf(t, a[i % Na], -v[i % Nv], 1 - w[i % Nw], sv[i % Nsv],
                          err[i % Nerr], switch_thresh, numf, sumf);
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
        if (t > 0 && isfinite(t)) { // sort response and calculate density
          if (out[i] == 1) { // response is "lower" so use unchanged parameters
              out[i] = denf(t, a[i % Na]/sigma[i % Nsig],
                            v[i % Nv]/sigma[i % Nsig], w[i % Nw],
                            sv[i % Nsv]/sigma[i % Nsig], err[i % Nerr],
                            switch_thresh, numf, sumf);
          } else { // response is "upper" so use alternate parameters
            out[i] = denf(t, a[i % Na]/sigma[i % Nsig],
                          -v[i % Nv]/sigma[i % Nsig], 1 - w[i % Nw],
                          sv[i % Nsv]/sigma[i % Nsig], err[i % Nerr],
                          switch_thresh, numf, sumf);
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
