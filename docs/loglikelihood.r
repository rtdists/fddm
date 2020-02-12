# Log-Likelihood functions for real data fitting





# Get initial parameter values
get_start_values <- function() {
  return(c(
    v1 = 0.1, # "upper" drift rate
    v0 = -0.1, # "lower" drift rate
    a = 1,
    w = 0.5
  ))
}




# Loglikelihood functions
ll_fs_eps_17 <- function(pars, rt, resp, truth, t0, eps) {
  rt1 <- rt[truth == 1]
  rt0 <- rt[truth == 0]
  resp1 <- resp[truth == 1]
  resp0 <- resp[truth == 0]

  # the truth is "upper" so use v1
  dens1 <- fs_eps_2017(rt = rt1, resp = resp1, v = pars[[1]],
                       a = pars[[3]], w = pars[[4]], t0 = t0, eps = eps)
  # the truth is "lower" so use v0
  dens0 <- fs_eps_2017(rt = rt0, resp = resp0, v = pars[[2]],
                       a = pars[[3]], w = pars[[4]], t0 = t0, eps = eps)

  densities <- c(dens1, dens0)
  if (any(densities == 0)) return(1e6)
  return(-sum(log(densities)))
}

ll_fs_BGK_17 <- function(pars, rt, resp, truth, t0, eps) {
  rt1 <- rt[truth == 1]
  rt0 <- rt[truth == 0]
  resp1 <- resp[truth == 1]
  resp0 <- resp[truth == 0]

  # the truth is "upper" so use v1
  dens1 <- fs_BGK_2017(rt = rt1, resp = resp1, v = pars[[1]],
                       a = pars[[3]], w = pars[[4]], t0 = t0, eps = eps)
  # the truth is "lower" so use v0
  dens0 <- fs_BGK_2017(rt = rt0, resp = resp0, v = pars[[2]],
                       a = pars[[3]], w = pars[[4]], t0 = t0, eps = eps)

  densities <- c(dens1, dens0)
  if (any(densities == 0)) return(1e6)
  return(-sum(log(densities)))
}

ll_fs_Nav_17 <- function(pars, rt, resp, truth, t0, eps) {
  rt1 <- rt[truth == 1]
  rt0 <- rt[truth == 0]
  resp1 <- resp[truth == 1]
  resp0 <- resp[truth == 0]

  # the truth is "upper" so use v1
  dens1 <- fs_Nav_2017(rt = rt1, resp = resp1, v = pars[[1]],
                       a = pars[[3]], w = pars[[4]], t0 = t0, eps = eps)
  # the truth is "lower" so use v0
  dens0 <- fs_Nav_2017(rt = rt0, resp = resp0, v = pars[[2]],
                       a = pars[[3]], w = pars[[4]], t0 = t0, eps = eps)

  densities <- c(dens1, dens0)
  if (any(densities == 0)) return(1e6)
  return(-sum(log(densities)))
}

ll_fl_Nav <- function(pars, rt, resp, truth, t0, eps) {
  rt1 <- rt[truth == 1]
  rt0 <- rt[truth == 0]
  resp1 <- resp[truth == 1]
  resp0 <- resp[truth == 0]

  # the truth is "upper" so use v1
  dens1 <- fl_Nav(rt = rt1, resp = resp1, v = pars[[1]],
                  a = pars[[3]], w = pars[[4]], t0 = t0, eps = eps)
  # the truth is "lower" so use v0
  dens0 <- fl_Nav(rt = rt0, resp = resp0, v = pars[[2]],
                  a = pars[[3]], w = pars[[4]], t0 = t0, eps = eps)

  densities <- c(dens1, dens0)
  if (any(densities == 0)) return(1e6)
  return(-sum(log(densities)))
}

ll_RWiener <- function(pars, rt, resp, truth, t0) {
  rt1 <- rt[truth == 1]
  rt0 <- rt[truth == 0]
  resp1 <- resp[truth == 1]
  resp0 <- resp[truth == 0]

  # the truth is "upper" so use v1
  dens1 <- dwiener(rt1, resp = resp1, delta = pars[[1]],
                   alpha = pars[[3]], beta = pars[[4]], tau = t0)
  # the truth is "lower" so use v0
  dens0 <- dwiener(rt0, resp = resp0, delta = pars[[2]],
                   alpha = pars[[3]], beta = pars[[4]], tau = t0)

  densities <- c(dens1, dens0)
  if (any(densities == 0)) return(1e6)
  return(-sum(log(densities)))
}

ll_RTDists <- function(pars, rt, resp, truth, t0) {
  rt1 <- rt[truth == 1]
  rt0 <- rt[truth == 0]
  resp1 <- resp[truth == 1]
  resp0 <- resp[truth == 0]

  # the truth is "upper" so use v1
  dens1 <- ddiffusion(rt1, resp1, v = pars[[1]],
                      a = pars[[3]], z = pars[[4]]*pars[[3]], t0 = t0)
  # the truth is "lower" so use v0
  dens0 <- ddiffusion(rt0, resp0, v = pars[[2]],
                      a = pars[[3]], z = pars[[4]]*pars[[3]], t0 = t0)

  densities <- c(dens1, dens0)
  if (any(densities == 0)) return(1e6)
  return(-sum(log(densities)))
}
