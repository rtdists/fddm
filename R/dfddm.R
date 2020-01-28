#' Density of Ratcliff Diffusion Decision Model
#' 
#' Density function for the Ratcliff diffusion decision model (DDM) with
#' following parameters: \code{a} (threshold separation), \code{w} (relative
#' starting point), \code{v} (drift rate), \code{t0} (non-decision time/response
#' time constant), and \code{sv} (inter-trial-variability of drift).
#' 
#' @param rt a vector of RTs. 
#' @param response ???
#' @param a threshold separation. Amount of information that is considered for a
#'   decision. Large values indicate a conservative decisional style. Typical
#'   range: 0.5 < \code{a} < 2
#' @param v drift rate. Average slope of the information accumulation process.
#'   The drift gives information about the speed and direction of the
#'   accumulation of information. Large (absolute) values of drift indicate a
#'   good performance. If received information supports the response linked to
#'   the upper threshold the sign will be positive and vice versa. Typical
#'   range: -5 < \code{v} < 5
#' @param t0 non-decision time or response time constant (in seconds). Lower
#'   bound for the duration of all non-decisional processes (encoding and
#'   response execution). Typical range: 0.1 < \code{t0} < 0.5
#' @param w relative starting point. Indicator of an a priori bias in decision
#'   making. When the relative starting point \code{w} deviates from
#'   \code{0.5}, the amount of information necessary for a decision differs
#'   between response alternatives. Default is \code{0.5} (i.e., no bias).
#' @param sv inter-trial-variability of drift rate. Standard deviation of a
#'   normal distribution with mean \code{v} describing the distribution of
#'   actual drift rates from specific trials. Values different from 0 can
#'   predict slow errors. Typical range: 0 < \code{sv} < 2. Default is 0.
#' @param log logical; if \code{TRUE}, probabilities p are given as log(p).
#' @param eps \code{numerical} scalar value. Precision of calculation.
#' @param n_terms_small Method for calculating number of terms used in the sum
#'   of the small time approximation of the DDM.
#' @param summation_small Method
#' @param scale some stuff
#' 
#' @details The different methods implemented
#' 
#' @references all the papers go here
#' 
#' @example examples/examples.diffusion.R
#' 
#' @return Density of the DDM with precision \code{eps} of length \code{rt}.
#' 
#' @useDynLib fddm, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @export
dfddm <- function(rt, response, 
                  a, v, t0, w = 0.5, 
                  sv = 0,
                  log = FALSE, 
                  n_terms_small = c("Foster", "Navarro", "Kesselmeier"),
                  summation_small = c("2017", "2014"),
                  scale = c("small", "large", "both"), 
                  eps = 0.000001
) {
  n_terms_small <- match.arg(n_terms_small)
  summation_small <- match.arg(summation_small)
  scale <- match.arg(scale)
  ## important bit here, maybe switch()
}
