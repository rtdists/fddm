// Second-Order Partial Derivatives of the Density function for the Ratcliff DDM

#include "helper_functions/fitting/partial_decs.h"
#include "helper_functions/fitting/declarations.h"




//' Second-Order Partial Derivative of 5-parameter DDM PDF with respect to v
//' (drift rate)
//'
//' Second-Order Partial Derivative of the density function for the 5-parameter
//'   variant of the Ratcliff diffusion decision model (DDM) with respect to v,
//'   the drift rate. This variant contains the following parameters:
//'   \code{v} (drift rate), \code{a} (threshold separation),
//'   \code{t0} (non-decision time/response time constant), \code{w}
//'   (relative starting point), \code{sv} (inter-trial variability of drift),
//'   and \code{sigma} (diffusion coefficient of underlying Wiener process).
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
//'   simply scales the parameters \code{v}, \code{a}, and \code{sv} as follows.
//'   See Details for more information. \itemize{ \item \code{a}
//'   \ifelse{html}{\out{&#8594;}}{\eqn{\to}} \code{a} \eqn{/} \code{sigma}
//'   \item \code{v} \ifelse{html}{\out{&#8594;}}{\eqn{\to}} \code{v} \eqn{/}
//'   \code{sigma} \item \code{sv} \ifelse{html}{\out{&#8594;}}{\eqn{\to}}
//'   \code{sv} \eqn{/} \code{sigma} }.
//'
//' @param sl_thresh Threshold for deciding when to use the "small-time" variant
//'   or the "large-time" variant. If the "effective response time" is greater
//'   than \code{sl_thresh} (i.e., \eqn{\frac{rt}{a^2} >} \code{sl_thresh}),
//'   then the "large-time" variant is used; otherwise, the "small-time" variant
//'   is used. Allowed values are any real number; however any non-positive
//'   number means that the "large-time" variant will always be used. Similarly,
//'   any very large positive number (e.g., +Inf) means that the "small-time"
//'   variant will always be used. Default value is \eqn{0.355}.
//'
//' @param err_tol Allowed error tolerance of the overall calculation. Since the
//'   partial derivative of the density function contains one infinite sum, this
//'   parameter defines the precision of the approximation to that infinite sum.
//'   If the provided error tolerance is less than \eqn{1e-300}, it is set to
//'   \eqn{1e-300}. Default is \eqn{1e-6}.
//'
//'
//'
//' @details
//'
//' All of the model inputs and parameters (\code{rt}, \code{response},
//' \code{v}, \code{a}, \code{t0}, \code{w}, \code{sv}, \code{sigma}) can be
//' input as a single value or as a vector of values. If input as a vector of
//' values, then the standard \code{R} input wrapping will occur.
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
//'
//'
//' @references
//'   Navarro, D. J., & Fuss, I. G. (2009). Fast and accurate
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
//'   Hartmann, R., Klauer, K. C. (2021). Partial derivatives for the
//'   first-passage time distribution in Wiener diffusion models. Journal of
//'   Mathematical Psychology, 103, 102550.
//'
//'   Ratcliff, R. (1978). A theory of memory retrieval. Psychological review,
//'   85(2), 59.
//'
//'
//'
//' @example examples/examples.pdf.partials.R
//'
//'
//'
//' @return A vector containing the second-order partial derivatives of the
//'   DDM PDF with precision \code{err_tol} whose length matches that of the
//'   longest input parameter (usually \code{rt}).
//'
//' @useDynLib fddm, .registration = TRUE
//' @import Rcpp
//' @export
// [[Rcpp::export]]
NumericVector dv2_dfddm(const NumericVector& rt,
                        const SEXP& response,
                        const NumericVector& v,
                        const NumericVector& a,
                        const NumericVector& t0,
                        const NumericVector& w = 0.5,
                        const NumericVector& sv = 0.0,
                        const NumericVector& sigma = 1.0,
                        const double& sl_thresh = 0.355, // no sum derivative
                        const NumericVector& err_tol = 0.000001)
{
  ParFunc parf = &dv2;
  return Rcpp::wrap(partial_pdf(parf, rt, response, v, a, t0, w, sv, sigma,
                                sl_thresh, err_tol));
}



//' Second-Order Partial Derivative of 5-parameter DDM PDF with respect to a
//' (threshold separation)
//'
//' Second-Order Partial Derivative of the density function for the 5-parameter
//'   variant of the Ratcliff diffusion decision model (DDM) with respect to a,
//'   the threshold separation. This variant contains the following parameters:
//'   \code{v} (drift rate), \code{a} (threshold separation),
//'   \code{t0} (non-decision time/response time constant), \code{w}
//'   (relative starting point), \code{sv} (inter-trial variability of drift),
//'   and \code{sigma} (diffusion coefficient of underlying Wiener process).
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
//'   simply scales the parameters \code{v}, \code{a}, and \code{sv} as follows.
//'   See Details for more information. \itemize{ \item \code{a}
//'   \ifelse{html}{\out{&#8594;}}{\eqn{\to}} \code{a} \eqn{/} \code{sigma}
//'   \item \code{v} \ifelse{html}{\out{&#8594;}}{\eqn{\to}} \code{v} \eqn{/}
//'   \code{sigma} \item \code{sv} \ifelse{html}{\out{&#8594;}}{\eqn{\to}}
//'   \code{sv} \eqn{/} \code{sigma} }.
//'
//' @param sl_thresh Threshold for deciding when to use the "small-time" variant
//'   or the "large-time" variant. If the "effective response time" is greater
//'   than \code{sl_thresh} (i.e., \eqn{\frac{rt}{a^2} >} \code{sl_thresh}),
//'   then the "large-time" variant is used; otherwise, the "small-time" variant
//'   is used. Allowed values are any real number; however any non-positive
//'   number means that the "large-time" variant will always be used. Similarly,
//'   any very large positive number (e.g., +Inf) means that the "small-time"
//'   variant will always be used. Default value is \eqn{0.5}.
//'
//' @param err_tol Allowed error tolerance of the overall calculation. Since the
//'   partial derivative of the density function contains the sum of two
//'   infinite sums, each approximation of these two infinite sums will have an
//'   individual error tolerance of \code{err_tol} / 2; thus the total overall
//'   error of the calculation will be at most \code{err_tol}. If the provided
//'   error tolerance is less than \eqn{1e-300}, it is set to \eqn{1e-300}.
//'   Default is \eqn{1e-6}.
//'
//'
//'
//' @details
//'
//' All of the model inputs and parameters (\code{rt}, \code{response},
//' \code{v}, \code{a}, \code{t0}, \code{w}, \code{sv}, \code{sigma}) can be
//' input as a single value or as a vector of values. If input as a vector of
//' values, then the standard \code{R} input wrapping will occur.
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
//'
//'
//' @references
//'   Navarro, D. J., & Fuss, I. G. (2009). Fast and accurate
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
//'   Hartmann, R., Klauer, K. C. (2021). Partial derivatives for the
//'   first-passage time distribution in Wiener diffusion models. Journal of
//'   Mathematical Psychology, 103, 102550.
//'
//'   Ratcliff, R. (1978). A theory of memory retrieval. Psychological review,
//'   85(2), 59.
//'
//'
//'
//' @example examples/examples.pdf.partials.R
//'
//'
//'
//' @return A vector containing the second-order partial derivatives of the
//'   DDM PDF with precision \code{err_tol} whose length matches that of the
//'   longest input parameter (usually \code{rt}).
//'
//' @useDynLib fddm, .registration = TRUE
//' @import Rcpp
//' @export
// [[Rcpp::export]]
NumericVector da2_dfddm(const NumericVector& rt,
                        const SEXP& response,
                        const NumericVector& v,
                        const NumericVector& a,
                        const NumericVector& t0,
                        const NumericVector& w = 0.5,
                        const NumericVector& sv = 0.0,
                        const NumericVector& sigma = 1.0,
                        const double& sl_thresh = 0.5, // with sum derivative
                        const NumericVector& err_tol = 0.000001)
{
  ParFunc parf = &da2;
  return Rcpp::wrap(partial_pdf(parf, rt, response, v, a, t0, w, sv, sigma,
                                sl_thresh, err_tol));
}



//' Second-Order Partial Derivative of 5-parameter DDM PDF with respect to t
//' (response time)
//'
//' Second-Order Partial Derivative of the density function for the 5-parameter
//'   variant of the Ratcliff diffusion decision model (DDM) with respect to t,
//'   the response time. This variant contains the following parameters:
//'   \code{v} (drift rate), \code{a} (threshold separation),
//'   \code{t0} (non-decision time/response time constant), \code{w}
//'   (relative starting point), \code{sv} (inter-trial variability of drift),
//'   and \code{sigma} (diffusion coefficient of underlying Wiener process).
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
//'   simply scales the parameters \code{v}, \code{a}, and \code{sv} as follows.
//'   See Details for more information. \itemize{ \item \code{a}
//'   \ifelse{html}{\out{&#8594;}}{\eqn{\to}} \code{a} \eqn{/} \code{sigma}
//'   \item \code{v} \ifelse{html}{\out{&#8594;}}{\eqn{\to}} \code{v} \eqn{/}
//'   \code{sigma} \item \code{sv} \ifelse{html}{\out{&#8594;}}{\eqn{\to}}
//'   \code{sv} \eqn{/} \code{sigma} }.
//'
//' @param sl_thresh Threshold for deciding when to use the "small-time" variant
//'   or the "large-time" variant. If the "effective response time" is greater
//'   than \code{sl_thresh} (i.e., \eqn{\frac{rt}{a^2} >} \code{sl_thresh}),
//'   then the "large-time" variant is used; otherwise, the "small-time" variant
//'   is used. Allowed values are any real number; however any non-positive
//'   number means that the "large-time" variant will always be used. Similarly,
//'   any very large positive number (e.g., +Inf) means that the "small-time"
//'   variant will always be used. Default value is \eqn{0.5}.
//'
//' @param err_tol Allowed error tolerance of the overall calculation. Since the
//'   partial derivative of the density function contains the sum of two
//'   infinite sums, each approximation of these two infinite sums will have an
//'   individual error tolerance of \code{err_tol} / 2; thus the total overall
//'   error of the calculation will be at most \code{err_tol}. If the provided
//'   error tolerance is less than \eqn{1e-300}, it is set to \eqn{1e-300}.
//'   Default is \eqn{1e-6}.
//'
//'
//'
//' @details
//'
//' All of the model inputs and parameters (\code{rt}, \code{response},
//' \code{v}, \code{a}, \code{t0}, \code{w}, \code{sv}, \code{sigma}) can be
//' input as a single value or as a vector of values. If input as a vector of
//' values, then the standard \code{R} input wrapping will occur.
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
//'
//'
//' @references
//'   Navarro, D. J., & Fuss, I. G. (2009). Fast and accurate
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
//'   Hartmann, R., Klauer, K. C. (2021). Partial derivatives for the
//'   first-passage time distribution in Wiener diffusion models. Journal of
//'   Mathematical Psychology, 103, 102550.
//'
//'   Ratcliff, R. (1978). A theory of memory retrieval. Psychological review,
//'   85(2), 59.
//'
//'
//'
//' @example examples/examples.pdf.partials.R
//'
//'
//'
//' @return A vector containing the second-order partial derivatives of the
//'   DDM PDF with precision \code{err_tol} whose length matches that of the
//'   longest input parameter (usually \code{rt}).
//'
//' @useDynLib fddm, .registration = TRUE
//' @import Rcpp
//' @export
// [[Rcpp::export]]
NumericVector dt2_dfddm(const NumericVector& rt,
                        const SEXP& response,
                        const NumericVector& v,
                        const NumericVector& a,
                        const NumericVector& t0,
                        const NumericVector& w = 0.5,
                        const NumericVector& sv = 0.0,
                        const NumericVector& sigma = 1.0,
                        const double& sl_thresh = 0.5, // with sum derivative
                        const NumericVector& err_tol = 0.000001)
{ // this is the same as dt02 because dt = -dt0, so dt^2 = dt0^2
  ParFunc parf = &dt02;
  return Rcpp::wrap(partial_pdf(parf, rt, response, v, a, t0, w, sv, sigma,
                                sl_thresh, err_tol));
}



//' Second-Order Partial Derivative of 5-parameter DDM PDF with respect to t0
//' (non-decision time)
//'
//' Second-Order Partial Derivative of the density function for the 5-parameter
//'   variant of the Ratcliff diffusion decision model (DDM) with respect to t0,
//'   the non-decision time. This variant contains the following parameters:
//'   \code{v} (drift rate), \code{a} (threshold separation),
//'   \code{t0} (non-decision time/response time constant), \code{w}
//'   (relative starting point), \code{sv} (inter-trial variability of drift),
//'   and \code{sigma} (diffusion coefficient of underlying Wiener process).
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
//'   simply scales the parameters \code{v}, \code{a}, and \code{sv} as follows.
//'   See Details for more information. \itemize{ \item \code{a}
//'   \ifelse{html}{\out{&#8594;}}{\eqn{\to}} \code{a} \eqn{/} \code{sigma}
//'   \item \code{v} \ifelse{html}{\out{&#8594;}}{\eqn{\to}} \code{v} \eqn{/}
//'   \code{sigma} \item \code{sv} \ifelse{html}{\out{&#8594;}}{\eqn{\to}}
//'   \code{sv} \eqn{/} \code{sigma} }.
//'
//' @param sl_thresh Threshold for deciding when to use the "small-time" variant
//'   or the "large-time" variant. If the "effective response time" is greater
//'   than \code{sl_thresh} (i.e., \eqn{\frac{rt}{a^2} >} \code{sl_thresh}),
//'   then the "large-time" variant is used; otherwise, the "small-time" variant
//'   is used. Allowed values are any real number; however any non-positive
//'   number means that the "large-time" variant will always be used. Similarly,
//'   any very large positive number (e.g., +Inf) means that the "small-time"
//'   variant will always be used. Default value is \eqn{0.5}.
//'
//' @param err_tol Allowed error tolerance of the overall calculation. Since the
//'   partial derivative of the density function contains the sum of two
//'   infinite sums, each approximation of these two infinite sums will have an
//'   individual error tolerance of \code{err_tol} / 2; thus the total overall
//'   error of the calculation will be at most \code{err_tol}. If the provided
//'   error tolerance is less than \eqn{1e-300}, it is set to \eqn{1e-300}.
//'   Default is \eqn{1e-6}.
//'
//'
//'
//' @details
//'
//' All of the model inputs and parameters (\code{rt}, \code{response},
//' \code{v}, \code{a}, \code{t0}, \code{w}, \code{sv}, \code{sigma}) can be
//' input as a single value or as a vector of values. If input as a vector of
//' values, then the standard \code{R} input wrapping will occur.
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
//'
//'
//' @references
//'   Navarro, D. J., & Fuss, I. G. (2009). Fast and accurate
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
//'   Hartmann, R., Klauer, K. C. (2021). Partial derivatives for the
//'   first-passage time distribution in Wiener diffusion models. Journal of
//'   Mathematical Psychology, 103, 102550.
//'
//'   Ratcliff, R. (1978). A theory of memory retrieval. Psychological review,
//'   85(2), 59.
//'
//'
//'
//' @example examples/examples.pdf.partials.R
//'
//'
//'
//' @return A vector containing the second-order partial derivatives of the
//'   DDM PDF with precision \code{err_tol} whose length matches that of the
//'   longest input parameter (usually \code{rt}).
//'
//' @useDynLib fddm, .registration = TRUE
//' @import Rcpp
//' @export
// [[Rcpp::export]]
NumericVector dt02_dfddm(const NumericVector& rt,
                         const SEXP& response,
                         const NumericVector& v,
                         const NumericVector& a,
                         const NumericVector& t0,
                         const NumericVector& w = 0.5,
                         const NumericVector& sv = 0.0,
                         const NumericVector& sigma = 1.0,
                         const double& sl_thresh = 0.5, // with sum derivative
                         const NumericVector& err_tol = 0.000001)
{ // this is the same as dt2 because dt = -dt0, so dt^2 = dt0^2
  ParFunc parf = &dt02;
  return Rcpp::wrap(partial_pdf(parf, rt, response, v, a, t0, w, sv, sigma,
                                sl_thresh, err_tol));
}



//' Second-Order Partial Derivative of 5-parameter DDM PDF with respect to w
//' (initial bias)
//'
//' Second-Order Partial Derivative of the density function for the 5-parameter
//'   variant of the Ratcliff diffusion decision model (DDM) with respect to w,
//'   the initial bias. This variant contains the following parameters:
//'   \code{v} (drift rate), \code{a} (threshold separation),
//'   \code{t0} (non-decision time/response time constant), \code{w}
//'   (relative starting point), \code{sv} (inter-trial variability of drift),
//'   and \code{sigma} (diffusion coefficient of underlying Wiener process).
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
//'   simply scales the parameters \code{v}, \code{a}, and \code{sv} as follows.
//'   See Details for more information. \itemize{ \item \code{a}
//'   \ifelse{html}{\out{&#8594;}}{\eqn{\to}} \code{a} \eqn{/} \code{sigma}
//'   \item \code{v} \ifelse{html}{\out{&#8594;}}{\eqn{\to}} \code{v} \eqn{/}
//'   \code{sigma} \item \code{sv} \ifelse{html}{\out{&#8594;}}{\eqn{\to}}
//'   \code{sv} \eqn{/} \code{sigma} }.
//'
//' @param sl_thresh Threshold for deciding when to use the "small-time" variant
//'   or the "large-time" variant in one half of the calculation of the partial
//'   derivative of the density function. This acts in exactly the same way as
//'   the parameter \code{max_terms_large} in the \code{dfddm()} function: it is
//'   the maximum number of terms to use for the "large-time" variant before
//'   switching to using the "small-time" variant. The other half of the
//'   calculation leverages a different type of comparison between the
//'   "large-time" and "small-time" variants that requires no input from the
//'   user. Allowed values are any real number; however any non-positive number
//'   means that the "small-time" variant will always be used. Similarly, any
//'   very large positive number (e.g., +Inf) means that the "large-time"
//'   variant will always be used. Default value is \eqn{1}.
//'
//' @param err_tol Allowed error tolerance of the overall calculation. Since the
//'   partial derivative of the density function contains the sum of two
//'   infinite sums, each approximation of these two infinite sums will have an
//'   individual error tolerance of \code{err_tol} / 2; thus the total overall
//'   error of the calculation will be at most \code{err_tol}. If the provided
//'   error tolerance is less than \eqn{1e-300}, it is set to \eqn{1e-300}.
//'   Default is \eqn{1e-6}.
//'
//'
//'
//' @details
//'
//' All of the model inputs and parameters (\code{rt}, \code{response},
//' \code{v}, \code{a}, \code{t0}, \code{w}, \code{sv}, \code{sigma}) can be
//' input as a single value or as a vector of values. If input as a vector of
//' values, then the standard \code{R} input wrapping will occur.
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
//'
//'
//' @references
//'   Navarro, D. J., & Fuss, I. G. (2009). Fast and accurate
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
//'   Hartmann, R., Klauer, K. C. (2021). Partial derivatives for the
//'   first-passage time distribution in Wiener diffusion models. Journal of
//'   Mathematical Psychology, 103, 102550.
//'
//'   Ratcliff, R. (1978). A theory of memory retrieval. Psychological review,
//'   85(2), 59.
//'
//'
//'
//' @example examples/examples.pdf.partials.R
//'
//'
//'
//' @return A vector containing the second-order partial derivatives of the
//'   DDM PDF with precision \code{err_tol} whose length matches that of the
//'   longest input parameter (usually \code{rt}).
//'
//' @useDynLib fddm, .registration = TRUE
//' @import Rcpp
//' @export
// [[Rcpp::export]]
NumericVector dw2_dfddm(const NumericVector& rt,
                        const SEXP& response,
                        const NumericVector& v,
                        const NumericVector& a,
                        const NumericVector& t0,
                        const NumericVector& w = 0.5,
                        const NumericVector& sv = 0.0,
                        const NumericVector& sigma = 1.0,
                        const double& sl_thresh = 1.0,
                        const NumericVector& err_tol = 0.000001)
{
  ParFunc parf = &dw2;
  return Rcpp::wrap(partial_pdf(parf, rt, response, v, a, t0, w, sv, sigma,
                                sl_thresh, err_tol));
}



//' Second-Order Partial Derivative of 5-parameter DDM PDF with respect to sv
//' (inter-trial variability in the drift rate)
//'
//' Second-Order Partial Derivative of the density function for the 5-parameter
//'   variant of the Ratcliff diffusion decision model (DDM) with respect to sv,
//'   the inter-trial variability in the drift rate. This variant contains the
//'   following parameters:
//'   \code{v} (drift rate), \code{a} (threshold separation),
//'   \code{t0} (non-decision time/response time constant), \code{w}
//'   (relative starting point), \code{sv} (inter-trial variability of drift),
//'   and \code{sigma} (diffusion coefficient of underlying Wiener process).
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
//'   simply scales the parameters \code{v}, \code{a}, and \code{sv} as follows.
//'   See Details for more information. \itemize{ \item \code{a}
//'   \ifelse{html}{\out{&#8594;}}{\eqn{\to}} \code{a} \eqn{/} \code{sigma}
//'   \item \code{v} \ifelse{html}{\out{&#8594;}}{\eqn{\to}} \code{v} \eqn{/}
//'   \code{sigma} \item \code{sv} \ifelse{html}{\out{&#8594;}}{\eqn{\to}}
//'   \code{sv} \eqn{/} \code{sigma} }.
//'
//' @param sl_thresh Threshold for deciding when to use the "small-time" variant
//'   or the "large-time" variant. If the "effective response time" is greater
//'   than \code{sl_thresh} (i.e., \eqn{\frac{rt}{a^2} >} \code{sl_thresh}),
//'   then the "large-time" variant is used; otherwise, the "small-time" variant
//'   is used. Allowed values are any real number; however any non-positive
//'   number means that the "large-time" variant will always be used. Similarly,
//'   any very large positive number (e.g., +Inf) means that the "small-time"
//'   variant will always be used. Default value is \eqn{0.355}.
//'
//' @param err_tol Allowed error tolerance of the overall calculation. Since the
//'   partial derivative of the density function contains one infinite sum, this
//'   parameter defines the precision of the approximation to that infinite sum.
//'   If the provided error tolerance is less than \eqn{1e-300}, it is set to
//'   \eqn{1e-300}. Default is \eqn{1e-6}.
//'
//'
//'
//' @details
//'
//' All of the model inputs and parameters (\code{rt}, \code{response},
//' \code{v}, \code{a}, \code{t0}, \code{w}, \code{sv}, \code{sigma}) can be
//' input as a single value or as a vector of values. If input as a vector of
//' values, then the standard \code{R} input wrapping will occur.
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
//'
//'
//' @references
//'   Navarro, D. J., & Fuss, I. G. (2009). Fast and accurate
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
//'   Hartmann, R., Klauer, K. C. (2021). Partial derivatives for the
//'   first-passage time distribution in Wiener diffusion models. Journal of
//'   Mathematical Psychology, 103, 102550.
//'
//'   Ratcliff, R. (1978). A theory of memory retrieval. Psychological review,
//'   85(2), 59.
//'
//'
//'
//' @example examples/examples.pdf.partials.R
//'
//'
//'
//' @return A vector containing the second-order partial derivatives of the
//'   DDM PDF with precision \code{err_tol} whose length matches that of the
//'   longest input parameter (usually \code{rt}).
//'
//' @useDynLib fddm, .registration = TRUE
//' @import Rcpp
//' @export
// [[Rcpp::export]]
NumericVector dsv2_dfddm(const NumericVector& rt,
                         const SEXP& response,
                         const NumericVector& v,
                         const NumericVector& a,
                         const NumericVector& t0,
                         const NumericVector& w = 0.5,
                         const NumericVector& sv = 0.0,
                         const NumericVector& sigma = 1.0,
                         const double& sl_thresh = 0.355, // no sum derivative
                         const NumericVector& err_tol = 0.000001)
{
  ParFunc parf = &dsv2;
  return Rcpp::wrap(partial_pdf(parf, rt, response, v, a, t0, w, sv, sigma,
                                sl_thresh, err_tol));
}
