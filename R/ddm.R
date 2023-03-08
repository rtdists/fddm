#' Estimation of 5-Parameter DDM
#'
#' Fit the 5-parameter DDM (Diffusion Decision Model) via maximum likelihood
#' estimation. The model for each DDM parameter can be specified symbolically
#' using R's formula interface. With the exception of the drift rate (which is
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
#' @details \code{ddm} uses \code{\link{model.matrix}} for transforming the
#'   symbolic description of the regression model underlying each parameter into
#'   estimated coefficients. The following provides a few examples:
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
#' @return Object of class \code{ddm} (i.e., a list with components as listed
#'   below) for which a number of common methods such as \code{print},
#'   \code{coef}, and \code{logLik} are implemented, see
#'   \code{\link{ddm-methods}}.
#'   \itemize{
#'     \item \code{coefficients} a named list whose elements are the values of
#'           the estimated model parameters
#'     \item \code{dpar} a character vector containing the names of the
#'           estimated model parameters
#'     \item \code{fixed_dpar} a named list whose elements are the values of the
#'           fixed model parameters
#'     \item \code{loglik} the value of the log-likelihood function at the
#'           optimized parameter values
#'     \item \code{hessians} a named list whose elements are the individual
#'           Hessians for each of the model parameters
#'     \item \code{vcov} a named list whose elements are the individual
#'           variance-covariance matrices for each of the model parameters
#'     \item \code{nobs} the number of observations in the data used for fitting
#'     \item \code{npar} the number of parameters used to fit the model (i.e.,
#'           the estimated model parameters plus any hyperparameters)
#'     \item \code{df.residual} the residual degrees of freedom (the number of
#'           observations - the number of parameters)
#'     \item \code{call} the original function call to \code{ddm}
#'     \item \code{formula} the formulas used in the model (\code{1} indicates
#'           that the model parameter was estimated with a single coefficient;
#'           \code{0} indicates that the model parameter was fixed)
#'     \item \code{dpar_formulas} a named list whose elements are the formulas
#'           for the model parameters (\code{1} indicates that the model
#'           parameter was estimated with a single coefficient; \code{0}
#'           indicates that the model parameter was fixed)
#'     \item \code{na.action} na.action
#'     \item \code{terms} a named list whose elements are , and whose last
#'           element is named \code{full} and shows the breakdown of the model
#'           with all model parameters
#'     \item \code{levels} a named list whose elements are the levels associated
#'           with any parameters that are factors (the elements are \code{NULL}
#'           if the parameter is not a factor), and whose last element is named
#'           \code{FULL} and shows all of the levels used in the model
#'     \item \code{contrasts} a named list whose elements are the type of
#'           contrasts used in the model
#'     \item \code{args_ddm} a named list whose elements are the optional
#'           arguments used in the calculation of the DDM log-likelihood
#'           function
#'     \item \code{link} a named list whose elements show information about the
#'           link function used for each model parameter (currently the only
#'           link function is the identity function)
#'     \item \code{converged} a logical indicating whether the optimization
#'           converged (\code{TRUE}) or not (\code{FALSE})
#'     \item \code{optim_info} a named list whose elements are information about
#'           the optimization process (e.g., the name of the algorithm used,
#'           the final value of the objective function, the number of
#'           evaluations of the gradient function, etc.)
#'     \item \code{model} the data used in the model (might need to check this)
#'     \item \code{response} the response data used in the model
#'     \item \code{mmatrix} a named list whose elements are the model matrices
#'           for each of the estimated parameters
#'     \item \code{compiled_model} C++ object that contains the compiled model
#'           (see list below for more details)
#'   }
#'   The C++ object accessible via the \code{compiled_model} component of the
#'   above R object of class \code{ddm} contains the following components:
#'   \itemize{
#'     \item \code{rt} a numeric vector of the response time data used in the
#'           model
#'     \item \code{response} an integer vector of the response data used in the
#'           model (coded such that \code{1} corresponds to the "lower" boundary
#'           and \code{2} corresponds to the "upper" boundary)
#'     \item \code{err_tol} the error tolerance used in the calculations for
#'           fitting the DDM
#'     \item \code{coefficients} a numeric vector containing the current set of
#'           coefficients for the formulas provided to the \code{ddm()} function
#'           call; the coefficients correspond to the DDM parameters in the
#'           following order: \code{v}, \code{a}, \code{t0}, \code{w}, \code{sv}
#'     \item \code{likelihood} a double containing the log-likelihood for the
#'           current set of \code{coefficients} (note this can be changed by
#'           calling the function \code{calculate_loglik()} below)
#'     \item \code{modmat_v} a numeric matrix containing the model matrix for
#'           \code{v}, the drift rate, determined by the formula input to the
#'           argument \code{drift} in the \code{ddm()} function call
#'     \item \code{modmat_a} a numeric matrix containing the model matrix for
#'           \code{a}, the boundary separation, determined by the formula input
#'           to the argument \code{boundary} in the \code{ddm()} function call
#'     \item \code{modmat_t0} a numeric matrix containing the model matrix for
#'           \code{t0}, the non-decision time, determined by the formula input
#'           to the argument \code{ndt} in the \code{ddm()} function call
#'     \item \code{modmat_w} a numeric matrix containing the model matrix for
#'           \code{w}, the inital bias, determined by the formula input to the
#'           argument \code{bias} in the \code{ddm()} function call
#'     \item \code{modmat_sv} a numeric matrix containing the model matrix for
#'           \code{sv}, the inter-trial variability in the drift rate,
#'           determined by the formula input to the argument \code{sv} in the
#'           \code{ddm()} function call
#'     \item \code{hess_v} a numeric matrix containing the Hessian for \code{v},
#'           the drift rate, whose dimensions are determined by the formula
#'           input to the argument \code{drift} in the \code{ddm()} function
#'           call
#'     \item \code{hess_a} a numeric matrix containing the Hessian for \code{a},
#'           the boundary separation, whose dimensions are determined by the
#'           formula input to the argument \code{drift} in the \code{ddm()}
#'           function call
#'     \item \code{hess_t0} a numeric matrix containing the Hessian for
#'           \code{t0}, the non-decision time, whose dimensions are determined
#'           by the formula input to the argument \code{drift} in the
#'           \code{ddm()} function call
#'     \item \code{hess_w} a numeric matrix containing the Hessian for \code{w},
#'           the initial bias, whose dimensions are determined by the formula
#'           input to the argument \code{drift} in the \code{ddm()} function
#'           call
#'     \item \code{hess_sv} a numeric matrix containing the Hessian for
#'           \code{sv}, the inter-trial variability in the drift rate, whose
#'           dimensions are determined by the formula input to the argument
#'           \code{drift} in the \code{ddm()} function call
#'     \item \code{vcov_v} a numeric matrix containing the variance-covariance
#'           matrix for \code{v}, the drift rate, whose dimensions are
#'           determined by the formula input to the argument \code{drift} in the
#'           \code{ddm()} function call
#'     \item \code{vcov_a} a numeric matrix containing the variance-covariance
#'           matrix for \code{a}, the boundary separation, whose dimensions are
#'           determined by the formula input to the argument \code{boundary} in
#'           the \code{ddm()} function call
#'     \item \code{vcov_t0} a numeric matrix containing the variance-covariance
#'           matrix for \code{t0}, the non-decision time, whose dimensions are
#'           determined by the formula input to the argument \code{ndt} in the
#'           \code{ddm()} function call
#'     \item \code{vcov_w} a numeric matrix containing the variance-covariance
#'           matrix for \code{w}, the inital bias, whose dimensions are
#'           determined by the formula input to the argument \code{bias} in the
#'           \code{ddm()} function call
#'     \item \code{vcov_sv} a numeric matrix containing the variance-covariance
#'           matrix for \code{v}, the inter-trial variability in the drift rate,
#'           whose dimensions are determined by the formula input to the
#'           argument \code{sv} in the \code{ddm()} function call
#'     \item \code{calculate_loglik} calculates and returns a double containing
#'           the negated log-likelihood (note that this will overwrite the
#'           \code{likelihood} component of the C++ object)
#'     \item \code{calculate_gradient} calculates and returns a numeric vector
#'           of the negated gradients for the provided coefficient values; the
#'           gradients are stored in the same manner as their corresponding
#'           \code{coefficents} (note that this will overwrite the
#'           \code{likelihood}) component of the C++ object)
#'     \item \code{calculate_hessians} calculates and returns a named list of
#'           the negated Hessians for each model parameter for the provided
#'           coefficient values (note that this will overwrite the
#'           \code{likelihood} component of the C++ object)
#'     \item \code{calculate_vcov} calculates and returns a named list of the
#'           variance-covariance matrices for each model parameter for the
#'           stored \code{coefficients}
#'     \item \code{calculate_standard_error} calculates and returns a numeric
#'           vector of the standard errors of the stored \code{coefficients};
#'           the standard errors are stored in the same manner as their
#'           corresponding \code{coefficients}
#'   }
#'
#' @importFrom stats .getXlevels make.link model.frame model.matrix nlminb terms
#' @importFrom methods new
#' @export
ddm <- function(drift, boundary = ~ 1, ndt = ~ 1, bias = 0.5, sv = 0,
                data,
                optim = "nlminb",
                args_optim = list(),
                args_ddm = list(err_tol = 1e-6), # update docs
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
    if (ncols[[par_name]] >= 1 && nrow(all_mm[[i]]) > 1) {
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
      #browser()
      if (par_name == "drift") {
        lm_mrt <- lm(rt ~ 0 + all_mm[[i]])
        #coef(lm_mrt)
        response_vec_num <- as.numeric(response_vec) - 1
        lm_acc <- lm(response_vec_num ~ 0 + all_mm[[i]])
        #coef(lm_acc)
        # tapply(data$response, INDEX = list(data$classification, data$difficulty), 
        #        FUN = function(x) mean(x == "blast"))
        var_rt <- var(rt)
        marg_par <- vector("numeric", ncols[[par_name]])
        new_start <- vector("numeric", ncols[[par_name]])
        for (j in seq_len(ncols[[par_name]])) {
          tmp_mm <- all_mm[[i]][ all_mm[[i]][,j] != 0, , drop = FALSE ]
          tmp_mrt <- sum(apply(tmp_mm[, seq_len(j), drop = FALSE], 2, mean) * coef(lm_mrt)[seq_len(j)])
          tmp_acc <- min(max(sum(apply(tmp_mm[, seq_len(j), drop = FALSE], 2, mean) * coef(lm_acc)[seq_len(j)]), 0), 1)
          
          # coef(lm_mrt)[seq_len(j)]
          # apply(tmp_mm[, seq_len(j), drop = FALSE], 2, mean)
          # ccc <- cumsum(apply(tmp_mm[, seq_len(j), drop = FALSE], 2, mean) * coef(lm_mrt)[seq_len(j)])
          # ccc[2] - ccc[1]
          # diff(rev(ccc * apply(tmp_mm[, seq_len(j), drop = FALSE], 2, mean)*-1))
          
          marg_par[j] <- ezddm(
            propCorrect = tmp_acc, 
            rtCorrectVariance_seconds = var_rt, 
            rtCorrectMean_seconds = tmp_mrt,
            s = 1, nTrials = nrow(tmp_mm)
          )[["v"]]
          
          if (j == 1) {
            new_start[j] <- marg_par[j]
          } else {
            new_start[j] <- diff(rev(-1*apply(tmp_mm[, seq_len(j), drop = FALSE], 2, mean) * marg_par[seq_len(j)]))
          }
          
          
        }
        browser()
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
  f <- new(fddm_fit, rt, response_vec, all_mm, args_ddm[["err_tol"]])

  #-------------------- Run Optimization --------------------------------------#
  fit_fun <- if (optim == "nlminb") fit_nlminb else optim
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

  ## prepare output object
  rval <- list(
    coefficients = opt$coefficients,
    dpar = names(formula_mm),
    fixed_dpar = all_ddm_constants,
    loglik = opt$loglik,
    hessians = f$calculate_hessians(all_coef)[names(all_ddm_formulas)],
    vcov = f$calculate_vcov()[names(all_ddm_formulas)],
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
