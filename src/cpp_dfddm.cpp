// Density function for the Ratcliff Diffusion Decision Model (DDM) PDF

#include "funcs.h"




// [[Rcpp::export]]
NumericVector cpp_dfddm(const NumericVector& rt,
                        const SEXP& response,
                        const NumericVector& a, const NumericVector& v,
                        const NumericVector& t0, const NumericVector& w,
                        const NumericVector& sv, const NumericVector& sigma,
                        const bool& log_prob,
                        const std::string& n_terms_small,
                        const std::string& summation_small,
                        const std::string& scale,
                        const int& max_terms_large,
                        const NumericVector& eps)
{
  // convert responses to false (0, lower) and true (1, upper)
  int Nres;
  vector<int> resp = convert_responses(response, Nres);



  // find Nmax (max length of parameter inputs)
  int Nrt  = rt.length();
  int Na   = a.length();
  int Nv   = v.length();
  int Nt0  = t0.length();
  int Nw   = w.length();
  int Nsv  = sv.length();
  int Nsig = sigma.length();
  int Neps = eps.length();
  int Nmax = max({Nrt, Nres, Na, Nv, Nt0, Nw, Nsv, Nsig, Neps});

  vector<double> a_c(Na);
  vector<double> t0_c(Nt0);
  vector<double> w_c(Nw);
  vector<double> sv_c(Nsv);
  vector<double> sigma_c(Nsig);
  vector<double> eps_c(Neps);

  // input checking
  if (!parameter_check(Nrt, Nres, Na, Nv, Nt0, Nw, Nsv, Nsig, Neps,
                       rt, a, t0, w, sv, sigma, eps,
                       a_c, t0_c, w_c, sv_c, sigma_c, eps_c)) {
    NumericVector empty_out(0);
    return empty_out;
  }



  // determine which method to use
  NumFunc numf;
  SumFunc sumf;
  DenFunc denf;
  double rt0;

  determine_method(n_terms_small, summation_small, scale,
                   numf, sumf, denf, rt0, log_prob);




  // loop through all inputs
  NumericVector out = calculate_pdf(Nrt, Nres, Na, Nv, Nt0, Nw, Nsv, Nsig, Neps,
                                    Nmax, rt, resp, a_c, v, t0_c, w_c, sv_c,
                                    sigma_c, eps_c, max_terms_large,
                                    numf, sumf, denf, rt0);




  return out;
}
