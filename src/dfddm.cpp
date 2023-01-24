// Density function for the Ratcliff Diffusion Decision Model (DDM)

#include "helper_functions/dfddm/declarations.h"



//' Density of Ratcliff Diffusion Decision Model
//'
//' Density function for the Ratcliff diffusion decision model (DDM) with
//' following parameters: \code{v} (drift rate), \code{a} (threshold
//' separation), \code{t0} (non-decision time/response time constant), \code{w}
//' (relative starting point), \code{sv} (inter-trial variability of drift), and
//' \code{sigma} (diffusion coefficient of underlying Wiener process). If you
//' are looking to fit the DDM, see [ddm()].
//'
//'
//'
//' @param rt A vector of response times (in seconds). If a response time is
//'   non-positve, then its density will evaluate to \eqn{0} if log = FALSE and
//'   \ifelse{html}{\out{-&#8734;}}{\eqn{-\infty}} if log = TRUE.
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
//' @param v Drift rate. Average slope of the information accumulation process.
//'   The drift gives information about the speed and direction of the
//'   accumulation of information. Large (absolute) values of drift indicate a
//'   good performance. If received information supports the response linked to
//'   the upper threshold, then the sign will be positive; similarly a negative
//'   value indicates that the received information supports the response linked
//'   to the lower threshold. Allowed range: \code{v} is a real number. Typical
//'   range: \eqn{-5 <} \code{v} \eqn{< 5}.
//'
//' @param a Threshold separation. Amount of information that is considered for
//'   a decision. Large values indicate a conservative decisional style. Allowed
//'   range: \eqn{0 <} \code{a}. Typical range: \eqn{0.5 <} \code{a} \eqn{< 2}.
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
//'   \code{sv} \eqn{/} \code{sigma} }.
//'
//' @param err_tol Allowed error tolerance of the density function. The density
//'   function contains an infinite sum that must be approximated; this
//'   parameter is the upper bound of the total error incurred from this
//'   approximation (in absolute value). See Details for more information.
//'   Default is \eqn{1e-6}.
//'
//' @param log Logical; if \code{TRUE}, probabilities \eqn{p} are given as
//'   \eqn{log(p)}. Default is \code{FALSE}, which gives the density on the
//'   probability scale.
//'
//' @param switch_mech Which switching mechanism to use in the choice of the
//'   "large-time" or "small-time" density function. Can be one of
//'   \{\code{"eff_rt"}, \code{"terms_large"}, \code{"terms"}, \code{"small"},
//'   \code{"large"}\}. Note that the large-time approximation is unstable for
//'   small effective response times (\eqn{(}\code{rt}\eqn{-}\code{t0}\eqn{)}
//'   \eqn{/(}\code{a}\eqn{*}\code{a}\eqn{) < 0.009}). See Details for more
//'   information. Default is \code{"eff_rt"}.
//'
//' @param switch_thresh Threshold for determining whether the effective
//'   response time (\eqn{(}\code{rt}\eqn{-}\code{t0}\eqn{)}
//'   \eqn{/(}\code{a}\eqn{*}\code{a}\eqn{)}) is "large" or "small". This
//'   parameter is only considered if \code{switch_mech = "eff_rt"} or
//'   \code{switch_mech = "terms_large"}.
//'   If \code{switch_mech = "eff_rt"}, an effective response time greater than
//'   \code{switch_thresh} is considered "large", and the "large-time" variant
//'   of the density function is used; otherwise, the "small-time" variant of
//'   the density function is used. The default is \eqn{0.8}.
//'   If \code{switch_mech = "terms_large"}, this parameter is treated as
//'   \eqn{ceiling(}\code{switch_thresh}\eqn{)}; the smallest integer that is
//'   not less than \code{switch_thresh}. In this case, the default is
//'   \eqn{ceiling(0.8) = 1}. See the \code{switch_mech} section of Details for
//'   more information.
//'   Note that if \code{switch_thresh}\eqn{ \le 0}, then the effective response
//'   time is always treated as "large"; contrarily, if \code{switch_thresh} =
//'   \ifelse{html}{\out{&#8734}}{\eqn{-\infty}} then the effective response
//'   time is always treated as "small". However, it is better to simply set
//'   \code{switch_mech = "large"} or \code{switch_mech = "small"} to always use
//'   the "large-time" or "small-time" variant, respectively.
//'
//' @param n_terms_small Which method to use for calculating the "small-time"
//'   approximation to the density function. Only applicable if
//'   \code{switch_mech = "terms"} or \code{switch_mech = "small"}; all other
//'   values of \code{switch_mech} cause this parameter to be ignored. If
//'   \code{switch_mech = "small"}, the allowed values are \code{"SWSE"},
//'   \code{"Gondan"}, and \code{"Navarro"}. The default value is \code{"SWSE"}.
//'   If \code{switch_mech = "terms"}, the allowed values are \code{"Gondan"}
//'   and \code{"Navarro"}. Note that if the user inputs
//'   \code{switch_mech = "terms"}, then the user must also explicitly input
//'   either \code{n_terms_small = "Gondan"} or
//'   \code{n_terms_small = "Navarro"}. See Details for more information.
//'
//' @param summation_small Which style of summation to use for the small-time
//'   approximation to the infinite sum. Can be one of \{\code{"2017"},
//'   \code{"2014"}\}. Only applicable if \code{switch_mech} is one of
//'   \{\code{"eff_rt"}, \code{"terms_large"}, \code{"terms"}, \code{"small"}\}.
//'   See Details for more information. Default is \code{"2017"}.
//'
//'
//'
//' @details
//'
//' All of the model inputs and parameters (\code{rt}, \code{response},
//' \code{v}, \code{a}, \code{t0}, \code{w}, \code{sv}, \code{sigma}) can be
//' input as a single value or as a vector of values. If input as a vector of
//' values, then the standard \code{R} recycling rules apply to ensure all 
//' inputs are of the same length.
//'
//' The default settings of \code{switch_mech = "eff_rt"},
//' \code{switch_thresh = "0.8"}, \code{n_terms_small = "SWSE"},
//' \code{summation_small = "2017"} produce the fastest and most accurate
//' results, as shown in our associated paper.
//'
//' \code{sv} - Both the "small-time" and "large-time" variants of the density
//' function have two further variants: one with a constant drift rate \code{v}
//' (i.e., \code{sv} \eqn{= 0}), and one with a variable drift rate \code{v}
//' (i.e., \code{sv} \eqn{> 0}). The details of the differences between these
//' two density functions can be found in our associated paper. To use the
//' density function with a constant drift rate, leave the parameter \code{sv}
//' to its default value of \code{sv = 0}, as this will indicate no drift to the
//' function. To use the density function with a variable drift rate, set the
//' parameter \code{sv} to some non-negative value, i.e., \code{sv} \eqn{> 0}.
//'
//' \code{sigma} - The default value of this parameter is \code{1} because it
//' only scales the parameters \code{v}, \code{a}, and \code{sv}, as shown
//' above. However, other formulations of the DDM may set \code{sigma = 0.1}
//' (see Ratcliff (1978), the fourth reference), so care must be taken when
//' comparing the results of different formulations.
//'
//' \code{err_tol} - The density function is composed of an infinite sum (that
//' must be approximated) and a multiplicative term outside the infinite sum,
//' \eqn{m}. The total error of the approximation is the error incurred from
//' truncating the infinite sum multiplied by \eqn{m}. Thus, to ensure that the
//' total error is bounded by the user-provided error tolerance, \code{err_tol},
//' the approximation to the infinite sum uses a modified error tolerance
//' (\code{err_tol}\eqn{ / m}). If the error tolerance is small and \eqn{m} is
//' large, the modified error tolerance can underflow to \eqn{0}. To protect
//' against this, we check that the modified error tolerance is at least
//' \eqn{1e-300} (near the smallest value that is representable by a floating
//' point number with double precision). If the modified error tolerance is
//' smaller than this threshold, then we silently change it to \eqn{1e-300}.
//' This case should only be encountered if extreme values of error tolerance
//' are used (i.e., on the scale of \eqn{1e-300}).
//'
//' \code{switch_mech} - The density function for the DDM has traditionally been
//' written in two forms: a "large-time" variant, and a "small-time" variant
//' (Navarro and Fuss, 2009). These two forms are more
//' efficient at calculating the density for large and small response times,
//' respectively. The parameter \code{switch_mech} determines how \code{dfddm}  
//' decides which of these two variants is used. 
//' \code{switch_mech = "small"} uses
//' only the "small-time" variant, and \code{switch_mech = "large"} uses only
//' the "large-time" variant. The "large-time" variant is unstable for small
//' effective response times (\eqn{(}\code{rt}\eqn{-}\code{t0}\eqn{)} \eqn{/
//' (}\code{a}\eqn{*}\code{a}\eqn{) < 0.009}) and may produce inaccurate
//' densities; thus, we do not recommend using the \code{switch_mech = "large"}
//' option if the inputs may contain such small effective response times. To
//' circumvent this accuracy issue and resolve the differing efficiencies of the
//' "large-time" and "small-time" variants, there are three switching mechanisms
//' that can be used to determine which of the two variants is more efficient.
//' First, \code{switch_mech = "terms"} is the traditional approach and
//' pre-calculates the number of terms required for the "large-time" and
//' "small-time" sums, and then uses whichever variant requires fewer terms.
//' This is the mechanism used in the Navarro and Fuss (2009) and Gondan,
//' Blurton, and Kesselmeier (2014) papers.
//' Second, \code{switch_mech = "terms_large"} pre-calculates the number of
//' terms only for the "large-time" variant, and compares that to the constant
//' value of \eqn{ceil(}\code{switch_thresh}\eqn{)} (default value of
//' \eqn{ceil(0.8) = 1}) to determine if the "large-time" variant is
//' sufficiently efficient.
//' Third, \code{switch_mech = "eff_rt"} determines whether a given effective
//' response time (\eqn{(}\code{rt}\eqn{-}\code{t0}\eqn{)}
//' \eqn{/(}\code{a}\eqn{*}\code{a}\eqn{)}) is considered "large" or "small" by
//' comparing it to the value of the parameter \code{switch_thresh} (default
//' value of \eqn{0.8}); it then uses the corresponding variant.
//' Both \code{switch_mech = "terms_large"} and \code{switch_mech = "eff_rt"}
//' only use the SWSE "small-time" approximation in the case that the
//' "small-time" variant is more efficient. \code{switch_mech = "eff_rt"} is the
//' most efficient method to determine which variant of the density function
//' should be used.
//'
//' \code{switch_thresh} - This parameter determines what effective response
//' times (\code{rt}\eqn{/(}\code{a} \eqn{*}\code{a}\eqn{)}) are "large" and
//' "small". The \code{paper_analysis} folder in the \code{fddm} GitHub
//' repository contains plots showing the relative efficiencies of a range of
//' values for the \code{switch_thresh} parameter when used with
//' \code{switch_mech = "eff_rt"} and also \code{switch_mech = "terms_large"}.
//' Note that this parameter changed name and purpose with the release of
//' \code{fddm} version 0.5-0.
//'
//' \code{n_terms_small} - The "small-time" variant has three different methods
//' for how to truncate the infinite sum in the density function. These
//' different methods are discussed extensively in our associated paper, but the
//' key distinction is that \code{n_terms_small = "SWSE"} uses a new method of
//' truncating the infinite sum. The \code{n_terms_small = "SWSE"} method is
//' currently recommended (when possible) because it is the fastest and most
//' stable algorithm when used with \code{switch_mech = "eff_rt"}.
//'
//' \code{summation_small} - The "large-time" variant of the density function
//' does not have any further variants, but the "small-time" variant has more
//' options with respect to evaluating the infinite sum. There are two
//' equivalent styles of summation, \code{summation_small = "2017"} and
//' \code{summation_small = "2014"}, of which the \code{"2017"} version
//' evaluates slightly faster and thus earns our recommendation. These different
//' styles of summation are discussed in our associated paper.
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
//' @example examples/examples.diffusion.R
//'
//'
//'
//' @return A vector containing the densities of the DDM with precision
//'   \code{err_tol} whose length matches that of the longest input parameter
//'   (usually \code{rt}).
//'
//' @seealso [ddm()]
//'
//' @useDynLib fddm, .registration = TRUE
//' @import Rcpp
//' @export
// [[Rcpp::export]]
NumericVector dfddm(const NumericVector& rt,
                    const SEXP& response,
                    const NumericVector& v,
                    const NumericVector& a,
                    const NumericVector& t0,
                    const NumericVector& w = 0.5,
                    const NumericVector& sv = 0.0,
                    const NumericVector& sigma = 1.0,
                    const NumericVector& err_tol = 0.000001,
                    const bool& log = false,
                    const std::string& switch_mech = "eff_rt",
                    double switch_thresh = 0.8,
                    const std::string& n_terms_small = "SWSE",
                    const std::string& summation_small = "2017")
{
  // determine which method to use
  NumFunc numf = NULL;
  SumFunc sumf = NULL;
  DenFunc denf = NULL;
  double rt0;

  determine_method(n_terms_small, summation_small, switch_mech,
                   switch_thresh, numf, sumf, denf, rt0, log);


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
  if (!parameter_check(Nrt, Nres, Nv, Na, Nt0, Nw, Nsv, Nsig, Nerr, Nmax,
                       rt, response, v, a, t0, w, sv, sigma, err_tol,
                       out, rt0)) {
    NumericVector empty_out(0);
    return empty_out;
  }


  // loop through all inputs, the vector `out` gets updated
  calculate_pdf(Nrt, Nv, Na, Nt0, Nw, Nsv, Nsig, Nerr,
                Nmax, rt, v, a, t0, w, sv, sigma,
                err_tol, out, switch_thresh,
                numf, sumf, denf, rt0);


  return Rcpp::wrap(out);
}
