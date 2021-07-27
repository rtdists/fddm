// Distribution function for the Ratcliff Diffusion Decision Model (DDM)

#include "parameter_checks.hpp"
#include "pfddm_helper_functions/func_defs.hpp"
#include "pfddm_helper_functions/helper_funcs.hpp"
#include "pfddm_helper_functions/distribution_funcs.hpp"
#include "pfddm_helper_functions/sum_funcs.hpp"
#include "pfddm_helper_functions/mills_funcs.hpp"



//' Distribution of Ratcliff Diffusion Decision Model
//'
//' Distribution function for the Ratcliff diffusion decision model (DDM) with
//' following parameters: \code{a} (threshold separation), \code{v} (drift
//' rate), \code{t0} (non-decision time/response time constant), \code{w}
//' (relative starting point), \code{sv} (inter-trial variability of drift), and
//' \code{sigma} (diffusion coefficient of underlying Wiener process).
//'
//'
//'
//' @param rt A vector of response times (in seconds). If a response time is
//'   non-positve, then its density will evaluate to \eqn{0} if log = FALSE and
//'   \ifelse{html}{\out{-<font style="vertical-align: middle;"
//'   size="5em">&#8734;</font>}}{\eqn{-\infty}} if log = TRUE.
//'
//' @param response Binary response(s) that correspond(s) to either the "lower"
//'   or "upper" threshold. This model parameter can either be a singular value
//'   or a vector. The value(s) in 'response' can be of the following data
//'   types:
//'     \itemize{
//'       \item integers or doubles (\eqn{1}
//'         \ifelse{html}{\out{&#8594;}}{\eqn{\to}} "lower", \eqn{2}
//'         \ifelse{html}{\out{&#8594;}}{\eqn{\to}} "upper");
//'       \item factors (the first level gets mapped to "lower", and the second
//'         level gets mapped to "upper"; any additional levels are ignored).
//'       \item strings (only the first character is checked, "L"
//'         \ifelse{html}{\out{&#8594;}}{\eqn{\to}} "lower" or "U"
//'         \ifelse{html}{\out{&#8594;}}{\eqn{\to}} "upper", case insensitive);
//'       \item logicals (FALSE \ifelse{html}{\out{&#8594;}}{\eqn{\to}} "lower",
//'         TRUE \ifelse{html}{\out{&#8594;}}{\eqn{\to}} "upper");
//'     }
//'
//' @param a Threshold separation. Amount of information that is considered for
//'   a decision. Large values indicate a conservative decisional style. Allowed
//'   range: \eqn{0 <} \code{a}. Typical range: \eqn{0.5 <} \code{a} \eqn{< 2}.
//'
//' @param v Drift rate. Average slope of the information accumulation process.
//'   The drift gives information about the speed and direction of the
//'   accumulation of information. Large (absolute) values of drift indicate a
//'   good performance. If received information supports the response linked to
//'   the upper threshold, then the sign will be positive; similarly a negative
//'   value indicates that the received information supports the response linked
//'   to the lower threshold. Allowed range: \code{v} is a real number. Typical
//'   range: \eqn{-5 <} \code{v} \eqn{< 5}.
//'
//' @param t0 Non-decision time or response time constant (in seconds). Lower
//'   bound for the duration of all non-decisional processes (encoding and
//'   response execution). If this value is greater than \code{rt}, then the
//'   resulting density is returned as if \code{rt} \eqn{ \le 0}. Allowed range:
//'   \eqn{0 \le} \code{t0}. Typical range: \eqn{0.1 <} \code{t0} \eqn{< 0.5}.
//'
//' @param w Relative starting point. Indicator of an a priori bias in decision
//'   making. When the relative starting point \code{w} deviates from \eqn{0.5},
//'   the amount of information necessary for a decision differs between
//'   response alternatives. Allowed range: \eqn{0 <} \code{w} \eqn{< 1}.
//'   Default value is \eqn{0.5} (i.e., no bias).
//'
//' @param sv Inter-trial-variability of drift rate. Standard deviation of a
//'   normal distribution with mean \code{v} describing the distribution of
//'   actual drift rates from specific trials. Values different from \eqn{0} can
//'   predict slow errors. Allowed range: \eqn{0 \le} \code{sv}. Typical range:
//'   \eqn{0 <} \code{sv} \eqn{< 2}. Default value is \eqn{0}, which indicates
//'   no drift in the function call. See Details for more information.
//'
//' @param sigma Diffusion coefficient of the underlying Wiener process. Allowed
//'   range: \eqn{0 <} \code{sigma}. Default value is \eqn{1}. This parameter
//'   simply scales the parameters \code{a}, \code{v}, and \code{sv} as follows.
//'   See Details for more information. \itemize{ \item \code{a}
//'   \ifelse{html}{\out{&#8594;}}{\eqn{\to}} \code{a} \eqn{/} \code{sigma}
//'   \item \code{v} \ifelse{html}{\out{&#8594;}}{\eqn{\to}} \code{v} \eqn{/}
//'   \code{sigma} \item \code{sv} \ifelse{html}{\out{&#8594;}}{\eqn{\to}}
//'   \code{sv} \eqn{/} \code{sigma} }
//'
//' @param log Logical; if \code{TRUE}, probabilities \eqn{p} are given as
//'   \eqn{log(p)}. Default is \code{FALSE}.
//'
//' @param err_tol Allowed error tolerance of the density function. Since the
//'   density function contains an infinite sum, this parameter defines the
//'   precision of the approximation to that infinite sum. Default is
//'   \eqn{1e-6}.
//'
//'
//'
//' @details
//'
//' All of the model inputs and parameters (\code{rt}, \code{response},
//' \code{a}, \code{v}, \code{t0}, \code{w}, \code{sv}, \code{sigma}) can be
//' input as a single value or as a vector of values. If input as a vector of
//' values, then the standard \code{R} input wrapping will occur.
//'
//' \code{sigma} - The default value of this parameter is \code{1} because it
//' only scales the parameters \code{a}, \code{v}, and \code{sv}, as shown
//' above. However, other formulations of the DDM may set \code{sigma = 0.1}
//' (see Ratcliff (1978), the fourth reference), so care must be taken when
//' comparing the results of different formulations.
//'
//'
//'
//' @references Navarro, D. J., & Fuss, I. G. (2009). Fast and accurate
//'   calculations for first-passage times in Wiener diffusion models. Journal
//'   of Mathematical Psychology, 53(4), 222-230.
//'
//'   Gondan, M., Blurton, S. P., & Kesselmeier, M. (2014). Even faster and even
//'   more accurate first-passage time densities and distributions for the
//'   Wiener diffusion model. Journal of Mathematical Psychology, 60, 20-22.
//'
//'   Blurton, S. P., Kesselmeier, M., & Gondan, M. (2017). The first-passage
//'   time distribution for the diffusion model with variable drift. Journal of
//'   Mathematical Psychology, 76, 7-12.
//'
//'   Ratcliff, R. (1978). A theory of memory retrieval. Psychological review,
//'   85(2), 59.
//'
//'
//'
//' @example examples/examples.distribution.R
//'
//'
//'
//' @return A vector containing the distribution of the DDM with precision
//'   \code{err_tol} whose length matches that of the longest input parameter
//'   (usually \code{rt}).
//'
//' @useDynLib fddm, .registration = TRUE
//' @import Rcpp
//' @export
// [[Rcpp::export]]
NumericVector pfddm(const NumericVector& rt,
                    const SEXP& response,
                    const NumericVector& a,
                    const NumericVector& v,
                    const NumericVector& t0,
                    const NumericVector& w = 0.5,
                    const NumericVector& sv = 0,
                    const NumericVector& sigma = 1,
                    const bool& log = 0,
                    const std::string& method = "1",
                    const NumericVector& err_tol = 0.000001)
{
  // determine which method to use (also log or non-log)
  DisFunc disf;
  double rt0;
  
  determine_method(method, disf, rt0, log);
  
  // determine lengths of parameter inputs, except response
  int Nrt  = rt.length();
  int Na   = a.length();
  int Nv   = v.length();
  int Nt0  = t0.length();
  int Nw   = w.length();
  int Nsv  = sv.length();
  int Nsig = sigma.length();
  int Nerr = err_tol.length();
  int Nmax = max({Nrt, Na, Nv, Nt0, Nw, Nsv, Nsig, Nerr}); // include Nres later
  int Nres;


  // initialize output, resized in convert_responses() inside parameter_check()
  vector<double> out;


  // check for invalid inputs, invalid inputs get marked in the vector `out`
  if (!parameter_check(Nrt, Nres, Na, Nv, Nt0, Nw, Nsv, Nsig, Nerr, Nmax,
                       rt, response, a, v, t0, w, sv, sigma, err_tol,
                       out, rt0)) {
    NumericVector empty_out(0);
    return empty_out;
  }

  // loop through all inputs, the vector `out` gets updated
  calculate_cdf(Nrt, Na, Nv, Nt0, Nw, Nsv, Nsig, Nerr,
                Nmax, rt, a, v, t0, w, sv, sigma,
                err_tol, out, rt0, disf);

  return Rcpp::wrap(out);
}
