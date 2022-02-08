// Helper functions for the density function dfddm



void ndetermine_method(const string& n_terms_small,
                      const string& summation_small,
                      const string& scale,
                      NumFunc& numf, SumFunc& sumf, DenFunc& denf,
                      double& rt0, const bool& log_prob)
{
  char n_terms_small0 = (!n_terms_small.empty()) ? n_terms_small[0] : EMPTYCHAR;
  char summation_small0 = (!summation_small.empty()) ?
    summation_small[summation_small.length()-1] : EMPTYCHAR;
  char scale0 = (!scale.empty()) ? scale[0] : EMPTYCHAR;

  if (log_prob) { // calculate log(probability)
    rt0 = -std::numeric_limits<double>::infinity();
    if (n_terms_small0 == 'S' || n_terms_small0 == 's') { // SWSE method
      if (scale0 == 'b' || scale0 == 'B') { // both
        denf = &nfc_log;
      } else if (scale0 == 's' || scale0 == 'S'){ // small
        denf = &nff_log;
      } else {
        stop("dfddm error: invalid function parameter 'scale': %s.", scale);
      }
      numf = NULL;
      if (summation_small0 == '7') { // 2017
        sumf = &nsmall_sum_eps_17;
      } else if (summation_small0 == '4') { // 2014
        sumf = &nsmall_sum_eps_14;
      } else {
        stop("dfddm error: invalid function parameter 'summation_small': %s.",
             summation_small);
      }
    } else {
      if (scale0 == 'l' || scale0 == 'L') { // large
        denf = &nfl_log;
        numf = NULL;
        sumf = NULL;
      } else {
        if (scale0 == 'b' || scale0 == 'B') { // both
          denf = &nfb_log;
        } else if (scale0 == 's' || scale0 == 'S') { // small
          denf = &nfs_log;
        } else {
          stop("dfddm error: invalid function parameter 'scale': %s.", scale);
        }
        if (n_terms_small0 == 'G' || n_terms_small0 == 'g') { // Gondan
          numf = &nks_Gon;
        } else if (n_terms_small0 == 'N' || n_terms_small0 == 'n') { // Navarro
          numf = &nks_Nav;
        } else {
          stop("dfddm error: invalid function parameter 'n_terms_small': %s.",
               n_terms_small);
        }
        if (summation_small0 == '7') { // 2017
          sumf = &nsmall_sum_2017;
        } else if (summation_small0 == '4') { // 2014
          sumf = &nsmall_sum_2014;
        } else {
          stop("dfddm error: invalid function parameter 'summation_small': %s.",
               summation_small);
        }
      }
    }
  } else { // calculate regular probability
    rt0 = 0;
    if (n_terms_small0 == 'S' || n_terms_small0 == 's') { // SWSE method
      if (scale0 == 'b' || scale0 == 'B') { // both
        denf = &nfc;
      } else if (scale0 == 's' || scale0 == 'S'){ // small
        denf = &nff;
      } else {
        stop("dfddm error: invalid function parameter 'scale': %s.", scale);
      }
      numf = NULL;
      if (summation_small0 == '7') { // 2017
        sumf = &nsmall_sum_eps_17;
      } else if (summation_small0 == '4') { // 2014
        sumf = &nsmall_sum_eps_14;
      } else {
        stop("dfddm error: invalid function parameter 'summation_small': %s.",
             summation_small);
      }
    } else {
      if (scale0 == 'l' || scale0 == 'L') { // large
        denf = &nfl;
        numf = NULL;
        sumf = NULL;
      } else {
        if (scale0 == 'b' || scale0 == 'B') { // both
          denf = &nfb;
        } else if (scale0 == 's' || scale0 == 'S') { // small
          denf = &nfs;
        } else {
          stop("dfddm error: invalid function parameter 'scale': %s.", scale);
        }
        if (n_terms_small0 == 'G' || n_terms_small0 == 'g') { // Gondan
          numf = &nks_Gon;
        } else if (n_terms_small0 == 'N' || n_terms_small0 == 'n') { // Navarro
          numf = &nks_Nav;
        } else {
          numf = &nks_Gon;
          stop("dfddm error: invalid function parameter 'n_terms_small': %s.",
               n_terms_small);
        }
        if (summation_small0 == '7') { // 2017
          sumf = &nsmall_sum_2017;
        } else if (summation_small0 == '4') { // 2014
          sumf = &nsmall_sum_2014;
        } else {
          sumf = &nsmall_sum_2017;
          stop("dfddm error: invalid function parameter 'summation_small': %s.",
               summation_small);
        }
      }
    }
  }
}



void ncalculate_pdf(const int& Nrt, const int& Na, const int& Nv, const int& Nt0,
                   const int& Nw, const int& Nsv, const int& Nsig,
                   const int& Nerr, const int& Nmax,
                   const NumericVector& rt,
                   const NumericVector& a, const NumericVector& v,
                   const NumericVector& t0, const NumericVector& w,
                   const NumericVector& sv, const NumericVector& sigma,
                   const NumericVector& err, vector<double>& out,
                   const double& sl_thresh,
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
                          err[i % Nerr], sl_thresh, numf, sumf);
          } else { // response is "upper" so use alternate parameters
            out[i] = denf(t, a[i % Na], -v[i % Nv], 1 - w[i % Nw], sv[i % Nsv],
                          err[i % Nerr], sl_thresh, numf, sumf);
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
                            sl_thresh, numf, sumf);
          } else { // response is "upper" so use alternate parameters
            out[i] = denf(t, a[i % Na]/sigma[i % Nsig],
                          -v[i % Nv]/sigma[i % Nsig], 1 - w[i % Nw],
                          sv[i % Nsv]/sigma[i % Nsig], err[i % Nerr],
                          sl_thresh, numf, sumf);
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
