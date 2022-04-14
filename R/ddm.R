#' Fitting function for the DDM
#'
#' This function fits the DDM to the data that is provided, using formulas for
#' the DDM parameters if desired.
#'
#' @param drift Formula for the drift rate (v)
#' @param boundary Formula for the boundary separation (a).
#' @param ntd Formula for the non-decision time (t0).
#' @param bias Formula for the relative initial bias (w).
#' @param sv Formula for the inter-trial variability in the drift rate (sv).
#'
#' @example examples/examples.ddm.R
#'
#' @return Currently, it's just the output from nlminb.
ddm <- function(drift, boundary = ~ 1, ndt = ~ 1, bias = 0.5, sv = 0,
                data,
                args_method = list(),
                args_ddm = list(err_tol = 1e-6, switch_thresh = 0.8),
                use_gradient = TRUE) {
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
  all_ddm_formulas <- all_ddm_pars[par_is_formula]
  all_vars <- unique(unlist(lapply(all_ddm_formulas, FUN = all.vars)))
  if (any(!all_vars %in% colnames(data))) {
    var_not_in_data <- all_vars[!all_vars %in% colnames(data)]
    stop("variables `", paste(var_not_in_data, collapse = ", "),
         "` not in data.", call. = FALSE)
  }

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
  # check NAs and remove (to prevent different Ns in different model matrices)
  data_use <- na.omit(data[, all_vars])
  if (nrow(data_use) != nrow(data)) {
    message("Removed ", nrow(data) - nrow(data_use), "observations with NAs, ",
            "N = ", nrow(data_use))
  }

  # extract RT and response and remove from formula
  rt_column <- all.vars(drift[[2]][[2]])
  rt <- data_use[[rt_column]]
  resp_column <- all.vars(drift[[2]][[3]])
  response <- data_use[[resp_column]]

  drift[[2]] <- NULL
  all_ddm_formulas[["drift"]] <- drift

  #-------------------- Create Model Matrices ---------------------------------#
  formula_mm <- lapply(all_ddm_formulas, model.matrix, data = data_use)
  miss_mm <- lapply(all_ddm_constants,
                    FUN = function(x) matrix(x, nrow = 1, ncol = 1))
  all_mm <- c(formula_mm, miss_mm)
  all_mm <- all_mm[names(all_ddm_pars)] # ensure correct order

  #-------------------- Check Estimability of Model Matrices ------------------#
  # also get initial values and bounds for optimization
  min_rt <- min(rt)
  inits <- list(
    "drift" = 0.0,
    "boundary" = 1.0,
    "ndt" = 0.33 * min_rt,
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
      par_rank <- qr(all_mm[[i]])[["rank"]]
      if (par_rank < ncols[[par_name]]) {
        all_mm[[i]] <- all_mm[[i]][, seq_len(par_rank), drop = FALSE]
        warning(paste("fddm_fit warning: model matrix for",
                      par_name, "was not estimable;",
                      ncols[[par_name]] - par_rank,
                      "columns were dropped from the right side.", sep = " "))
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

  #-------------------- Create fddm_fit Object --------------------------------#
  f <- new(fddm_fit, rt, response, all_mm,
           args_ddm[["err_tol"]], args_ddm[["switch_thresh"]])

  #-------------------- Run Optimization --------------------------------------#
  if (use_gradient) {
    store_this <- nlminb(init_vals, objective = f$calculate_loglik,
       gradient = f$calculate_gradient,
       lower = lower_bds, upper = upper_bds,
       control = args_method
       )
  } else {
    store_this <- nlminb(init_vals, objective = f$calculate_loglik,
       lower = lower_bds, upper = upper_bds,
       control = args_method
       )
  }

  return(store_this)
}
