
fit_nlminb <- function(init,
                       objective, gradient,
                       lower, upper,
                       control,
                       ...) ## ... arguments are ignored
{
  out <- nlminb(start = init,
                objective = objective,
                gradient = gradient,
                lower = lower, upper = upper,
                control = control)
  list(
    coefficients = out$par,
    loglik = -out$objective,
    converged = out$convergence == 0,
    optim = out
  )
}

fit_optim <- function(init,
                      objective, gradient,
                      lower, upper,
                      method,
                      control,
                      ...) ## ... arguments are ignored
{
  out <- optim(par = init,
               fn = objective,
               gr = gradient,
               lower = lower, upper = upper,
               method = method,
               control = control)
  list(
    coefficients = out$par,
    loglik = -out$value,
    converged = out$convergence == 0,
    optim = out
  )
}

ez_init_vals <- function(ncoeffs, min_rt, formula_mm, ncols, rt, response_vec) {
  init_vals <- numeric(ncoeffs)
  init_pars <- c(
    "drift" = 0.0,
    "boundary" = 1.0,
    "ndt" = 0.5 * min_rt,
    "bias" = 0.5,
    "sv" = 0.1,
    "diff" = 0.0
  )
  ez_parnames <- c(
    "drift" = "v",
    "boundary" = "a",
    "ndt" = "Ter"
  )

  for (i in seq_along(formula_mm)) { # only doing this for not fixed params
    parname <- names(formula_mm[i])
    idx <- sum(ncols[seq_len(i - 1)])
    if (parname %in% names(ez_parnames)) { # drift, boundary, ndt
      tmp_iv <- numeric(ncols[i])
      for (j in seq_len(ncols[i])) {
        tmp_mm <- formula_mm[[i]][formula_mm[[i]][, j] != 0, , drop = FALSE]
        weights <- apply(tmp_mm[, seq_len(j), drop = FALSE], 2, mean)
        if (weights[j] == 0) { # coef is a dif with zero weight
          init_vals[idx + j] <- init_pars[["diff"]]
        } else { # coef could be an int or a diff with nonzero weight
          # get initial values using EZ-Diffusion (Wagenmakers et al. 2007)
          lm_mrt <- lm(rt ~ 0 + formula_mm[[i]])
          lm_acc <- lm((as.numeric(response_vec) - 1) ~ 0 + formula_mm[[i]])
          var_rt <- var(rt)
          tmp_mrt <- sum(weights * coef(lm_mrt)[seq_len(j)])
          tmp_acc <- min(max(sum(weights * coef(lm_acc)[seq_len(j)]), 0), 1)
          marg_par_j <- ezddm(
            propCorrect = tmp_acc,
            rtCorrectVariance_seconds = var_rt,
            rtCorrectMean_seconds = tmp_mrt,
            s = 1, nTrials = nrow(tmp_mm)
          )[[ez_parnames[[parname]]]]
          tmp_iv[j] <- (marg_par_j - sum(weights[-j] * tmp_iv[seq_len(j-1)])
                        ) / weights[j]
          if (weights[j] == 1 && all(weights[-j] == 0)) { # coef is an int
            if (parname == "ndt") { # check for ezddm weirdness
                if (tmp_iv[j] < 0) {
                  warning("ddm warning: the initial value for ndt is ",
                  "negative, so it will be changed to 0.1 * min(rt) = ",
                  0.1 * min_rt)
                  tmp_iv[j] <- 0.1 * min_rt
                } else if (tmp_iv[j] >= min_rt) {
                  warning("ddm warning: the initial value for ndt is ",
                  "greater than min(rt), so it will be changed to ",
                  "0.9 * min(rt) = ", 0.9 * min_rt, " < min(rt) = ", min_rt)
                  tmp_iv[j] <- 0.9 * min_rt
                }
              }
          }
        }
      }
      init_vals[(idx + 1):(idx + j)] <- tmp_iv
    } else { # bias, sv
      for (j in seq_len(ncols[i])) {
        tmp_mm <- formula_mm[[i]][formula_mm[[i]][, j] != 0, , drop = FALSE]
        weights <- apply(tmp_mm[, seq_len(j), drop = FALSE], 2, mean)
        if (weights[j] == 1 && all(weights[-j] == 0)) { # coef is an int
          init_vals[idx + j] <- init_pars[[parname]]
        } else { # coef is a dif
          init_vals[idx + j] <- init_pars[["diff"]]
        }
      }
    }
  }
  return(init_vals)
}

random_init_vals <- function(ncoeffs, min_rt, formula_mm, ncols) {
  init_vals <- numeric(ncoeffs)
  init_pars <- c(
    "drift" = rnorm(n = 1, mean = 0, sd = 1.5),
    "boundary" = rgamma(n = 1, shape = 20, rate = 12),
    "ndt" = min_rt * rbeta(n = 1, shape1 = 30, shape2 = 30),
    "bias" = rbeta(n = 1, shape1 = 30, shape2 = 30),
    "sv" = rgamma(n = 1, shape = 2, rate = 4),
    "diff" = rnorm(n = 1, mean = 0, sd = 0.05)
  )

  for (i in seq_along(formula_mm)) { # only doing this for not fixed params
    parname <- names(formula_mm[i])
    idx <- sum(ncols[seq_len(i - 1)])
    for (j in seq_len(ncols[i])) {
      tmp_mm <- formula_mm[[i]][formula_mm[[i]][, j] != 0, , drop = FALSE]
      weights <- apply(tmp_mm[, seq_len(j), drop = FALSE], 2, mean)
      if (weights[j] == 1 && all(weights[-j] == 0)) { # coef is an int
        init_vals[idx + j] <- init_pars[[parname]]
      } else { # coef is a dif
        init_vals[idx + j] <- init_pars[["diff"]]
      }
    }
  }

  return(init_vals)
}
