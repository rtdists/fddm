#' Density of Ratcliff Diffusion Decision Model
#'
#' Density function for the Ratcliff diffusion decision model (DDM) with
#' following parameters: \code{a} (threshold separation), \code{v} (drift rate),
#' \code{t0} (non-decision time/response time constant), \code{w} (relative
#' starting point), and \code{sv} (inter-trial variability of drift).
#'
#'
#'
#' @param rt A vector of response times (in seconds). If a response time is
#'   non-positve, then its density will evaluate to \eqn{0} if log = FALSE and
#'   \eqn{-Inf} if log = TRUE.
#'
#' @param response Binary response(s) that correspond(s) to either the "lower"
#'   or "upper" threshold. This model parameter can either be a singular value
#'   or a vector. The value(s) in 'response' can be of the following data types:
#'   \itemize{ \item logicals (FALSE \ifelse{html}{\out{&#8594;}}{\eqn{\to}}
#'   "lower", TRUE \ifelse{html}{\out{&#8594;}}{\eqn{\to}} "upper"); \item
#'   integers or doubles (\eqn{1} \ifelse{html}{\out{&#8594;}}{\eqn{\to}}
#'   "lower", \eqn{2} \ifelse{html}{\out{&#8594;}}{\eqn{\to}} "upper"); \item
#'   strings (only the first character is checked, "L"
#'   \ifelse{html}{\out{&#8594;}}{\eqn{\to}} "lower" or "U"
#'   \ifelse{html}{\out{&#8594;}}{\eqn{\to}} "upper", case insensitive); \item
#'   factors (the first level gets mapped to "lower", and the second level gets
#'   mapped to "upper"; any additional levels are ignored). }
#'
#' @param a Threshold separation. Amount of information that is considered for a
#'   decision. Large values indicate a conservative decisional style. Allowed
#'   range: \eqn{0 <} \code{a}. Typical range: \eqn{0.5 <} \code{a} \eqn{< 2}.
#'
#' @param v Drift rate. Average slope of the information accumulation process.
#'   The drift gives information about the speed and direction of the
#'   accumulation of information. Large (absolute) values of drift indicate a
#'   good performance. If received information supports the response linked to
#'   the upper threshold, then the sign will be positive; similarly a negative
#'   value indicates that the received information supports the response linked
#'   to the lower threshold. Allowed range: \code{v} is a real number. Typical
#'   range: \eqn{-5 <} \code{v} \eqn{< 5}.
#'
#' @param t0 Non-decision time or response time constant (in seconds). Lower
#'   bound for the duration of all non-decisional processes (encoding and
#'   response execution). If this value is greater than \code{rt}, then the
#'   resulting density is returned as if \code{rt} \eqn{ \leq 0}. Allowed range:
#'   \eqn{0 <} \code{t0}. Typical range: \eqn{0.1 <} \code{t0} \eqn{< 0.5}.
#'
#' @param w Relative starting point. Indicator of an a priori bias in decision
#'   making. When the relative starting point \code{w} deviates from \eqn{0.5},
#'   the amount of information necessary for a decision differs between response
#'   alternatives. Allowed range: \eqn{0 <} \code{w} \eqn{< 1}. Default is
#'   \eqn{0.5} (i.e., no bias).
#'
#' @param sv Inter-trial-variability of drift rate. Standard deviation of a
#'   normal distribution with mean \code{v} describing the distribution of
#'   actual drift rates from specific trials. Values different from \eqn{0} can
#'   predict slow errors. Allowed range: \eqn{0 <} \code{sv}. Typical range:
#'   \eqn{0 <} \code{sv} \eqn{< 2}. Default is \eqn{-1}, which indicates no
#'   drift in the function call. See Details for more information.
#'
#' @param log Logical; if \code{TRUE}, probabilities \eqn{p} are given as
#'   \eqn{log(p)}. Default is \code{FALSE}.
#'
#' @param n_terms_small Which method for calculating number of terms used in the
#'   approximation of the infinite sum in the small-time approximation of the
#'   density function. Can be one of \{\code{"SWSE"}, \code{"Gondan"},
#'   \code{"Navarro"}\}. Only applicable if \code{scale} is one of
#'   \{\code{"small"}, \code{"both"}\}. See Details for more information.
#'   Default is \code{"Gondan"}.
#'
#' @param summation_small Which style of summation to use for the small-time
#'   approximation to the infinite sum. Can be one of \{\code{"2017"},
#'   \code{"2014"}\}. Only applicable if \code{scale} is one of
#'   \{\code{"small"}, \code{"both"}\}. See Details for more information.
#'   Default is \code{"2017"}.
#'
#' @param scale Which density function to use. Can be one of \{\code{"small"},
#'   \code{"large"}, \code{"both"}\}. Note that the large-time approximation is
#'   unstable for small effective response times (\code{rt}\eqn{/(}\code{a}
#'   \eqn{*}\code{a}\eqn{) < 0.009}). See Details for more information. Default
#'   is \code{"both"}.
#'
#' @param err_tol Allowed error tolerance of the density function. Since the
#'   density function contains an infinite sum, this parameter defines the
#'   precision of the approximation to that infinite sum. Default is \eqn{1e-6}.
#'
#'
#'
#' @details
#'
#'   All of the model inputs and parameters (\code{rt}, response, \code{a},
#'   \code{v}, \code{t0}, \code{w}, and \code{sv}) can be input as a single value
#'   or as a vector of values. If input as a vector of values, then the standard
#'   \code{R} input wrapping will occur.
#'
#'   The default settings of \code{n_terms_small = "Gondan"},
#'   code{summation_small = "2017"}, \code{scale = "both"} produce the fastest
#'   and most accurate results, as shown in our associated paper.
#'
#'   \code{scale} - The density function for the DDM has traditionally been
#'   written in two forms: a "large-time" variant, and a "small-time" variant.
#'   The parameter \code{scale} determines which of these variants will be used
#'   in the calculation; \code{scale = "large"} uses the "large-time" variant,
#'   and \code{scale = "small"} uses the "small-time" variant. The "large-time"
#'   variant is unstable for small effective response times
#'   (\code{rt}\eqn{/(}\code{a}\eqn{*}\code{a}\eqn{) < 0.009}) and produces
#'   inaccurate densities; thus we do not recommend using only the \code{scale =
#'   "large"} option if the inputs contain such small response times. To
#'   circumvent this issue, the \code{scale = "both"} option utilizes both the
#'   "small-time" and "large-time" variants by determining which variant is more
#'   computationally efficient before calculating the density. Even though the
#'   "large-time" density function is often significantly slower than the
#'   "small-time" variant, it is extremely efficient in some areas of the
#'   parameter space, and so the \code{scale = "both"} option is the fastest.
#'
#'   \code{sv} - Both the "small-time" and "large-time" variants of the density
#'   function have two further variants: one with a constant drift rate \code{v}
#'   (i.e., \code{sv} \eqn{= 0}), and one with a variable drift rate \code{v}
#'   (i.e., \code{sv} \eqn{> 0}). The details of the differences between these
#'   two density functions can be found in our associated paper. To use the
#'   density function with a constant drift rate, leave the parameter \code{sv}
#'   to its default value of \code{sv = -1}, as this will indicate no drift to
#'   the function. To use the density function with a variable drift rate, set
#'   the parameter \code{sv} to some positive value, e.g. \code{sv = .75}.
#'
#'   \code{summation_small} - The "large-time" variant of the density function
#'   does not have any further variants, but the "small-time" has more options
#'   with respect to evaluating the infinite sum. There are two equivalent
#'   styles of summation, \code{summation_small = "2017"} and
#'   \code{summation_small = "2014"}, of which the \code{"2017"} version
#'   evaluates slightly faster and thus earns our recommendation. These
#'   different styles of summation are discussed in our associated paper.
#'
#'   \code{n_terms_small} - The "small-time" variant also has three different
#'   methods for how to truncate the infinite sum in the density function. These
#'   different methods are discussed extensively in our associated paper, but
#'   the key distinction is that \code{n_terms_small = "Foster"} uses a new
#'   method of truncating the infinite sum. The \code{n_terms_small = "Gondan"}
#'   method is currently recommended because it is the most stable algorithm
#'   when used with \code{scale = "both"}.
#'
#'
#'
#' @references Navarro, D. J., & Fuss, I. G. (2009). Fast and accurate
#'   calculations for first-passage times in Wiener diffusion models. Journal of
#'   Mathematical Psychology, 53(4), 222-230.
#'
#'   Gondan, M., Blurton, S. P., & Kesselmeier, M. (2014). Even faster and even
#'   more accurate first-passage time densities and distributions for the Wiener
#'   diffusion model. Journal of Mathematical Psychology, 60, 20-22.
#'
#'   Blurton, S. P., Kesselmeier, M., & Gondan, M. (2017). The first-passage
#'   time distribution for the diffusion model with variable drift. Journal of
#'   Mathematical Psychology, 76, 7-12.
#'
#'
#'
#' @example examples/examples.diffusion.R
#'
#'
#'
#' @return A vector containing the densities of the DDM with precision
#'   \code{err_tol} whose length matches that of the longest input parameter
#'   (usually \code{rt}).
#'
#' @useDynLib fddm, .registration = TRUE
#' @export
dfddm <- function(rt, response,
                  a, v, t0, w = 0.5,
                  sv = -1,
                  log = FALSE,
                  n_terms_small = "Gondan",
                  summation_small = "2017",
                  scale = "both",
                  err_tol = 0.000001)
{
  return(cpp_dfddm(rt, response, a, v, t0, w, sv, log,
                   n_terms_small, summation_small, scale, err_tol))
}
