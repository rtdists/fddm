#' Estimation of 5-Parameter DDM
#'
#' Fit the 5-parameter DDM (Diffusion Decision Model) via maximum likelihood
#' estimation. The model for each DDM parameter can be specified symbolically
#' using R's formulas interface. With the exception of the drift rate (which is
#' always estimated) all parameters can be either fixed or estimates.
#'
#' @param drift Two-sided formula. The left-hand side describes the response,
#'   the right-hand side provides a symbolic description of the regression model
#'   underlying the drift rate (v). The left-hand side needs to specify the
#'   column in the data containing response time and corresponding binary
#'   decision, concatenated by \code{+} , e.g., \code{rt + response ~ ...}
#' @param boundary,ndt,bias,sv Either a one-sided formula providing a symbolic
#'   description of the regression model or a scalar number given the value this
#'   parameter should be fixed to. Boundary separation (a), non-decision time
#'   (t0), relative initial bias (w), or inter-trial variability in the drift
#'   rate (sv)
#' @param data,na.action,subset arguments controlling formula processing via
#'   \code{\link{model.frame}}.
#' @param optim character string or fitting function indicating which numerical
#'   optimisation method should be used. The default \code{"nlminb"} uses the
#'   corresponding function.
#' @param args_optim additional control arguments passed to \code{control}
#'   argument of optimisation function specified in \code{optim}.
#' @param args_ddm additional arguments passed to density function.
#' @param use_gradient logical. Should gradient be used during numerical
#'   optimisation? Default is \code{TRUE}.
#' @param compiled_model,model,mmatrix,response logicals. If \code{TRUE} the
#'   corresponding components of the fit (the compiled model object, the model
#'   frame, the model matrix, the response matrix) are returned.
#' @param contrasts optional list. See the contrasts.arg of
#'   \code{\link{model.matrix.default}}
#' 
#' @details \code{ddm} uses \code{\link{model.matrix}} for transforming the symbolic description of the regression model underlying each parameter into estimated coefficients. The following provides a few examples:
#' 
#' \itemize{
#'   \item \code{~ 1} estimates a single coefficient, named \code{(Intercept)}
#'   \item \code{~ condition} estimates the intercept plus k - 1 coefficients
#'   for a factor with k levels (e.g., intercept plus one coefficient if
#'   condition has two levels). The interpretation of the coefficients depend on
#'   the factor contrasts employed, which are usually based on the contrasts
#'   specified in \code{options("contrasts")}. For the default \code{treatment}
#'   contrasts (\code{\link{contr.treatment}}), the intercept corresponds to the
#'   first factor level and the additional coefficients correspond to the
#'   difference from the intercept (i.e., first factor level). When using
#'   \code{contr.sum} the intercept correspond to the grand mean and the
#'   additional coefficients correspond to the differences from the grand mean.
#'   \item code{~ 0 + condition} estimates no intercept but one coefficient per
#'   factor level. This specification can also be used to get one coefficient
#'   per cell for a multi-factorial design, e.g., \code{~ 0 +
#'   condition1:condition2}.
#'   \item code{~ 0 + condition1 + condition1:condition2} estimates one
#'   "intercept" per level of \code{condition1} factor (which is not called
#'   intercept) plus k - 1 difference parameters from the condition-specific
#'   intercept for the k-levels of \code{condition2}. The interpretation of the
#'   difference parameters again depends on the contrasts used (e.g., treatment
#'   vs. sum-to-zero contrasts, see examples). This formula specification can
#'   often make sense for the drift rate when \code{condition1} is the factor
#'   (such as item type) mapped to upper and lower response boundary of the DDM
#'   and \code{condition2} is another factor by which we want the drift rate to
#'   differ. This essentially gives one overall drift rate per response boundary
#'   plus the differences from this overall one (note again that with treatment
#'   contrasts this overall drift rate is the first factor level of
#'   \code{condition2}).
#' }
#' 
#' To get meaningful results it is necessary to estimate separate drift rates
#' for the different condition/item-types that are mapped onto the upper and
#' lower boundary of the diffusion model.
#' 
#' If a non-default fitting function is used, it needs to minimise the negative
#' log-likelihood, accept the following arguments, \code{init, objective,
#' gradient, lower, upper, control} , and return a list with the following
#' arguments \code{coefficients, loglik, converged, optim} (where
#' \code{converged} is boolean and \code{optim} can be an arbitrary list with
#' additional information).
#'
#' @example examples/examples.ddm.R
#'
#' @return Object of class \code{ddm} for which a number of common methods such
#'   as \code{print}, \code{coef}, and \code{logLik} are implemented, see
#'   \code{\link{ddm-methods}}.
#'
#' @importFrom stats .getXlevels make.link model.frame model.matrix nlminb terms
#' @importFrom methods new
#' @export
ddm <- function(drift, boundary = ~ 1, ndt = ~ 1, bias = 0.5, sv = 0,
                data,
                optim = "nlminb",
                args_optim = list(),
                args_ddm = list(err_tol = 1e-6, switch_thresh = 0.8),
                use_gradient = TRUE,
                compiled_model = TRUE, model = TRUE,
                mmatrix = TRUE, response = TRUE,
                na.action, subset,
                contrasts = NULL) {
  cl <- match.call()

  ## prepare single formula for model.frame

  #-------------------- Extract and Check Formulas ----------------------------#
  all_ddm_pars <- list(
    drift = drift,
    boundary = boundary,
    ndt = ndt,
    bias = bias,
    sv = sv
  )
  # is parameter formula or not?
  par_is_formula <- vapply(all_ddm_pars, FUN = inherits, FUN.VALUE = NA,
                           what = "formula")
  # check formulas:
  all_ddm_formulas <- all_ddm_pars
  all_ddm_formulas[!par_is_formula] <- 
    lapply(all_ddm_formulas[!par_is_formula], 
           FUN = function(x) ~ 0)
  # all_ddm_formulas <- all_ddm_pars[par_is_formula]

  ## uses Formula package: https://cran.r-project.org/package=Formula
  full_formula <- do.call(Formula::as.Formula, args = unname(all_ddm_formulas))
  # attr(full_formula, "ddm_parameters") <- names(all_ddm_formulas)

  ## see: https://developer.r-project.org/model-fitting-functions.html
  mf <- match.call(expand.dots = FALSE)
  mf$formula <- full_formula
  m <- match(c("formula", "data", "subset", "weights", "na.action",
                 "offset"), names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- quote(stats::model.frame) 
  mf <- eval.parent(mf)
  mt <- attr(mf, "terms")

  # check fixed parameters:
  all_ddm_constants <- all_ddm_pars[!par_is_formula]
  length_constants <- vapply(all_ddm_constants,
                             FUN = length, FUN.VALUE = NA_integer_)
  if (any(length_constants != 1)) {
    stop("Fixed parameter `", which(length_constants != 1),
         "` must be of length 1.", call. = FALSE)
  }
  # more checks? maybe later

  #-------------------- Extract and Check Data --------------------------------#
  # extract RT and response and remove from formula
  response_df <- Formula::model.part(
    object = full_formula,
    data = mf,
    lhs = 1,
    drop = TRUE
  )

  rt_column <- all.vars(drift[[2]][[2]])
  rt <- response_df[[rt_column]]
  resp_column <- all.vars(drift[[2]][[3]])
  response_vec <- response_df[[resp_column]]

  drift[[2]] <- NULL
  mm_formulas <- all_ddm_formulas[par_is_formula]
  mm_formulas[["drift"]] <- drift

  #-------------------- Create Model Matrices ---------------------------------#
  formula_mm <- vector("list", length(full_formula)[2])
  names(formula_mm) <- names(all_ddm_pars)
  rterms <- vector("list", length(full_formula)[2])
  names(rterms) <- names(all_ddm_pars)
  for (i in seq_along(formula_mm)) {
    rterms[[i]] <- terms(full_formula, lhs = 0, rhs = i)
    formula_mm[[i]] <- model.matrix(
      object = rterms[[i]],
      data = mf,
      contrasts.arg = contrasts
    )
  }
  all_mm <- formula_mm
  all_mm[!par_is_formula] <- lapply(all_ddm_constants,
                    FUN = function(x) matrix(x, nrow = 1, ncol = 1))

  ### link function (currently only uses identity, so no argument yet)
  all_links <- list(
    drift = "identity",
    boundary = "identity",
    ndt = "identity",
    bias = "identity",
    sv = "identity"
  )
  links <- lapply(all_links[names(all_ddm_formulas)], make.link)

  #-------------------- Check Estimability of Model Matrices ------------------#
  # also get initial values and bounds for optimization
  min_rt <- min(rt)
  inits <- list(
    "drift" = 0.0,
    "boundary" = 1.0,
    "ndt" = 0.9 * min_rt,
    "bias" = 0.5,
    "sv" = 0.0,
    "zero" = 0.0
  )
  l_bds <- list(
    "drift" = -Inf,
    "boundary" = 0.0,
    "ndt" = 0.0,
    "bias" = 0.0,
    "sv" = 0.0,
    "diff" = -Inf
  )
  u_bds <- list(
    "drift" = Inf,
    "boundary" = Inf,
    "ndt" = min_rt,
    "bias" = 1.0,
    "sv" = Inf,
    "diff" = Inf
  )

  ncols <- lapply(all_mm, FUN = ncol)
  init_vals <- numeric()
  lower_bds <- numeric()
  upper_bds <- numeric()
  for (i in seq_along(all_mm)) {
    par_name <- names(all_mm[i])[1]
    if (ncols[[par_name]] >= 1 & nrow(all_mm[[i]]) > 1) {
      # check estimability
      rank_warn <- FALSE
      # par_rank <- qr(all_mm[[i]])[["rank"]]
      while (qr(all_mm[[i]])[["rank"]] < ncol(all_mm[[i]])) {
        all_mm[[i]] <- all_mm[[i]][, seq_len(qr(all_mm[[i]])[["rank"]]),
                                   drop = FALSE]
        rank_warn <- TRUE
      }
      if (rank_warn) {
        warning("model matrix for ", par_name, " was rank deficient; ",
                ncols[[par_name]] - ncol(all_mm[[i]]),
                " column(s) dropped from right side.",
                call. = FALSE)
        ncols[[par_name]] <- ncol(all_mm[[i]])
      }
      # get initial values and bounds for estimation
      if ("(Intercept)" == colnames(all_mm[[i]])[1]) {
        init_vals <- c(init_vals, inits[[par_name]],
                      rep(inits[["zero"]], ncols[[par_name]] - 1))
        lower_bds <- c(lower_bds, l_bds[[par_name]],
                      rep(l_bds[["diff"]], ncols[[par_name]] - 1))
        upper_bds <- c(upper_bds, u_bds[[par_name]],
                      rep(u_bds[["diff"]], ncols[[par_name]] - 1))
      } else {
        init_vals <- c(init_vals, rep(inits[[par_name]], ncols[[par_name]]))
        lower_bds <- c(lower_bds, rep(l_bds[[par_name]], ncols[[par_name]]))
        upper_bds <- c(upper_bds, rep(u_bds[[par_name]], ncols[[par_name]]))
      }
    }
  }
  formula_mm <- all_mm[par_is_formula]

  #-------------------- Create fddm_fit Object --------------------------------#
  f <- new(fddm_fit, rt, response_vec, all_mm,
           args_ddm[["err_tol"]], args_ddm[["switch_thresh"]])

  #-------------------- Run Optimization --------------------------------------#
  if (optim == "nlminb") fit_fun <- fit_nlminb
  else fit_fun <- optim
  opt <- fit_fun(
    init = init_vals,
    objective = f$calculate_loglik,
    gradient = if (use_gradient) f$calculate_gradient else NULL,
    lower = lower_bds, upper = upper_bds,
    control = args_optim)

  #-------------------- Prepare Output ----------------------------------------#
  ## prepare coefficient list:
  all_coef <- opt$coefficients
  names(all_coef) <- unlist(lapply(formula_mm,  colnames))
  coef_lengths <- vapply(formula_mm,  ncol, 0)
  coef_map <- rep(names(coef_lengths), coef_lengths)
  coef_map <- factor(coef_map, levels = unique(coef_map))
  opt$coefficients <- split(all_coef, coef_map)

  coef_list <- vector("list", length(all_ddm_formulas))
  names(coef_list) <- names(all_ddm_formulas)

  ## calculate Hessian and variance-covariance matrix
  f$calculate_hessians(all_coef)
  hess_temp <- list(
    drift = f$hess_v,
    boundary = f$hess_a,
    ndt = f$hess_t0,
    bias = f$hess_w,
    sv = f$hess_sv)
  f$calculate_vcov()
  vcov_temp <- list(
    drift = f$vcov_v,
    boundary = f$vcov_a,
    ndt = f$vcov_t0,
    bias = f$vcov_w,
    sv = f$vcov_sv)

  ## prepare output object
  rval <- list(
    coefficients = opt$coefficients,
    dpar = names(formula_mm),
    fixed_dpar = all_ddm_constants,
    loglik = opt$loglik,
    hessians = hess_temp[names(all_ddm_formulas)],
    vcov = vcov_temp[names(all_ddm_formulas)],
    nobs = nrow(response_df),
    npar = length(init_vals),
    df.residual = nrow(response_df) - length(init_vals)
  )
  rval$call <- cl
  rval$formula <- full_formula
  rval$dpar_formulas <- all_ddm_formulas
  rval$na.action <- attr(mf, "na.action")
  rval$terms <- c(rterms, full = mt)
  rval$levels <- c(lapply(rterms, .getXlevels, m = mf),
                   full = list(.getXlevels(mt, mf)))
  rval$contrasts <- lapply(formula_mm, attr, which = "contrasts")
  rval$args_ddm <- args_ddm
  rval$link <- links
  rval$converged <- opt$converged
  rval$optim_info <- list(
    optim = optim,
    args_optim = args_optim,
    use_gradient = use_gradient,
    value = opt$optim
  )
  if (compiled_model) {
    rval$compiled_model <- f
  }
  if (model) {
    rval$model <- mf
  }
  if (response) {
    rval$response <- Formula::model.part(full_formula, mf,
                                         lhs = 1, rhs = 0)
  }
  if (mmatrix) {
    rval$mmatrix <- formula_mm
  }
  class(rval) <- "ddm"
  return(rval)
}
