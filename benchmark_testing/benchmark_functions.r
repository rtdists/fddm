library("microbenchmark")
library("rtdists")
library("RWiener")
source("inst/extdata/Kesselmeier_density.R")



####################### Constant Drift Rate
bm <- function(RT, resp, V, A, W, t0 = 1e-4, err_tol = 1e-6,
               times = 1000, unit = "us")
{
  if (resp < 0.5) {
    rresp <- "lower"
  } else {
    rresp <- "upper"
  }

  nf <- 14 # number of functions being benchmarked
  nRT <- length(RT) # number of response times
  nV <- length(V) # number of drift rates
  nA <- length(A) # number of boundary separations
  nW <- length(W) # number of starting points
  fnames <- c("fs_Fos_17", "fs_Fos_14", "fs_Kes_17", "fs_Kes_14",
              "fs_Nav_17", "fs_Nav_14", "fl_Nav_09",
              "fb_Kes_17", "fb_Kes_14", "fb_Nav_17", "fb_Nav_14",
              "rtdists", "RWiener", "Kesselmeier")

  # Initialize the dataframe to contain the microbenchmark results
  mbm_res <- data.frame(matrix(ncol=times+5, nrow=nf*nRT*nV*nA*nW))
  colnames(mbm_res) <- c('RT', 'V', 'A', 'W', 'FuncName',
                         paste0("bm", seq_len(times)))
  start <- 1
  stop <- nf
  sl_time <- seq_len(times)+5

  # Loop through each combination of parameters and record microbenchmark results
  for (rt in 1:nRT) {
    for (v in 1:nV) {
      for (a in 1:nA) {
        for (w in 1:nW) {
          mbm <- microbenchmark(
            fs_Fos_17 = dfddm(rt = RT[rt], response = resp, a = A[a],
                              v = V[v], t0 = t0, w = W[w],
                              log = FALSE, n_terms_small = "Foster",
                              summation_small = "2017", scale = "small",
                              err_tol = err_tol),
            fs_Fos_14 = dfddm(rt = RT[rt], response = resp, a = A[a],
                              v = V[v], t0 = t0, w = W[w],
                              log = FALSE, n_terms_small = "Foster",
                              summation_small = "2014", scale = "small",
                              err_tol = err_tol),
            fs_Kes_17 = dfddm(rt = RT[rt], response = resp, a = A[a],
                              v = V[v], t0 = t0, w = W[w],
                              log = FALSE, n_terms_small = "Kesselmeier",
                              summation_small = "2017", scale = "small",
                              err_tol = err_tol),
            fs_Kes_14 = dfddm(rt = RT[rt], response = resp, a = A[a],
                              v = V[v], t0 = t0, w = W[w],
                              log = FALSE, n_terms_small = "Kesselmeier",
                              summation_small = "2014", scale = "small",
                              err_tol = err_tol),
            fs_Nav_17 = dfddm(rt = RT[rt], response = resp, a = A[a],
                              v = V[v], t0 = t0, w = W[w],
                              log = FALSE, n_terms_small = "Navarro",
                              summation_small = "2017", scale = "small",
                              err_tol = err_tol),
            fs_Nav_14 = dfddm(rt = RT[rt], response = resp, a = A[a],
                              v = V[v], t0 = t0, w = W[w],
                              log = FALSE, n_terms_small = "Navarro",
                              summation_small = "2014", scale = "small",
                              err_tol = err_tol),
            fl_Nav_09 = dfddm(rt = RT[rt], response = resp, a = A[a],
                              v = V[v], t0 = t0, w = W[w],
                              log = FALSE, n_terms_small = NULL,
                              summation_small = NULL, scale = "large",
                              err_tol = err_tol),
            fb_Kes_17 = dfddm(rt = RT[rt], response = resp, a = A[a],
                              v = V[v], t0 = t0, w = W[w],
                              log = FALSE, n_terms_small = "Kesselmeier",
                              summation_small = "2017", scale = "both",
                              err_tol = err_tol),
            fb_Kes_14 = dfddm(rt = RT[rt], response = resp, a = A[a],
                              v = V[v], t0 = t0, w = W[w],
                              log = FALSE, n_terms_small = "Kesselmeier",
                              summation_small = "2014", scale = "both",
                              err_tol = err_tol),
            fb_Nav_17 = dfddm(rt = RT[rt], response = resp, a = A[a],
                              v = V[v], t0 = t0, w = W[w],
                              log = FALSE, n_terms_small = "Navarro",
                              summation_small = "2017", scale = "both",
                              err_tol = err_tol),
            fb_Nav_14 = dfddm(rt = RT[rt], response = resp, a = A[a],
                              v = V[v], t0 = t0, w = W[w],
                              log = FALSE, n_terms_small = "Navarro",
                              summation_small = "2014", scale = "both",
                              err_tol = err_tol),
            rtdists = ddiffusion(RT[rt], rresp, a = A[a], v = V[v],
                                 t0 = t0, z = W[w]*A[a]),
            RWiener = dwiener(RT[rt], resp = rresp, alpha = A[a],
                              delta = V[v], tau = t0, beta = W[w],
                              give_log = FALSE),
            Kesselmeier = fs14_R(t = RT[rt]-t0, a = A[a], v = V[v],
                                 w = W[w], eps = err_tol), # only "lower" resp
            times = times, unit = unit)
          # add the rt, v, a, w values and function names to the dataframe
          mbm_res[start:stop, 1] <- rep(RT[rt], nf)
          mbm_res[start:stop, 2] <- rep(V[v]  , nf)
          mbm_res[start:stop, 3] <- rep(A[a]  , nf)
          mbm_res[start:stop, 4] <- rep(W[w]  , nf)
          mbm_res[start:stop, 5] <- fnames
          # add the microbenchmark results to the dataframe
          for (i in 1:nf) {
            bm <- subset(mbm, expr==fnames[i])
            mbm_res[start+i-1, sl_time] <- bm$time
          }
          # iterate start and stop values
          start = start + nf
          stop = stop + nf
        }
      }
    }
  }
  return(mbm_res)
}


bm_vec <- function(RT, resp, V, A, W, t0 = 1e-4, err_tol = 1e-6,
                   times = 1000, unit = "us")
{
  if (resp < 0.5) {
    rresp <- "lower"
  } else {
    rresp <- "upper"
  }

  nf <- 14 # number of functions being benchmarked
  nV <- length(V) # number of drift rates
  nA <- length(A) # number of boundary separations
  nW <- length(W) # number of starting points
  rresp <- rep(rresp, length(RT)) # for RWiener
  fnames <- c("fs_Fos_17", "fs_Fos_14", "fs_Kes_17", "fs_Kes_14",
              "fs_Nav_17", "fs_Nav_14", "fl_Nav_09",
              "fb_Kes_17", "fb_Kes_14", "fb_Nav_17", "fb_Nav_14",
              "rtdists", "RWiener", "Kesselmeier")

  # Initialize the dataframe to contain the microbenchmark results
  mbm_res <- data.frame(matrix(ncol=times+4, nrow=nf*nV*nA*nW))
  colnames(mbm_res) <- c('V', 'A', 'W', 'FuncName',
                         paste0("bm", seq_len(times)))
  start <- 1
  stop <- nf
  sl_time <- seq_len(times)+4

  # Loop through each combination of parameters and record microbenchmark results
  for (v in 1:nV) {
    for (a in 1:nA) {
      for (w in 1:nW) {
        mbm <- microbenchmark(
          fs_Fos_17 = dfddm(rt = RT, response = resp, a = A[a],
                            v = V[v], t0 = t0, w = W[w],
                            log = FALSE, n_terms_small = "Foster",
                            summation_small = "2017", scale = "small",
                            err_tol = err_tol),
          fs_Fos_14 = dfddm(rt = RT, response = resp, a = A[a],
                            v = V[v], t0 = t0, w = W[w],
                            log = FALSE, n_terms_small = "Foster",
                            summation_small = "2014", scale = "small",
                            err_tol = err_tol),
          fs_Kes_17 = dfddm(rt = RT, response = resp, a = A[a],
                            v = V[v], t0 = t0, w = W[w],
                            log = FALSE, n_terms_small = "Kesselmeier",
                            summation_small = "2017", scale = "small",
                            err_tol = err_tol),
          fs_Kes_14 = dfddm(rt = RT, response = resp, a = A[a],
                            v = V[v], t0 = t0, w = W[w],
                            log = FALSE, n_terms_small = "Kesselmeier",
                            summation_small = "2014", scale = "small",
                            err_tol = err_tol),
          fs_Nav_17 = dfddm(rt = RT, response = resp, a = A[a],
                            v = V[v], t0 = t0, w = W[w],
                            log = FALSE, n_terms_small = "Navarro",
                            summation_small = "2017", scale = "small",
                            err_tol = err_tol),
          fs_Nav_14 = dfddm(rt = RT, response = resp, a = A[a],
                            v = V[v], t0 = t0, w = W[w],
                            log = FALSE, n_terms_small = "Navarro",
                            summation_small = "2014", scale = "small",
                            err_tol = err_tol),
          fl_Nav_09 = dfddm(rt = RT, response = resp, a = A[a],
                            v = V[v], t0 = t0, w = W[w],
                            log = FALSE, n_terms_small = NULL,
                            summation_small = NULL, scale = "large",
                            err_tol = err_tol),
          fb_Kes_17 = dfddm(rt = RT, response = resp, a = A[a],
                            v = V[v], t0 = t0, w = W[w],
                            log = FALSE, n_terms_small = "Kesselmeier",
                            summation_small = "2017", scale = "both",
                            err_tol = err_tol),
          fb_Kes_14 = dfddm(rt = RT, response = resp, a = A[a],
                            v = V[v], t0 = t0, w = W[w],
                            log = FALSE, n_terms_small = "Kesselmeier",
                            summation_small = "2014", scale = "both",
                            err_tol = err_tol),
          fb_Nav_17 = dfddm(rt = RT, response = resp, a = A[a],
                            v = V[v], t0 = t0, w = W[w],
                            log = FALSE, n_terms_small = "Navarro",
                            summation_small = "2017", scale = "both",
                            err_tol = err_tol),
          fb_Nav_14 = dfddm(rt = RT, response = resp, a = A[a],
                            v = V[v], t0 = t0, w = W[w],
                            log = FALSE, n_terms_small = "Navarro",
                            summation_small = "2014", scale = "both",
                            err_tol = err_tol),
          rtdists = ddiffusion(RT, rresp, a = A[a], v = V[v],
                               t0 = t0, z = W[w]*A[a]),
          RWiener = dwiener(RT, resp = rresp, alpha = A[a],
                            delta = V[v], tau = t0, beta = W[w],
                            give_log = FALSE),
          Kesselmeier = fs14_R(t = RT-t0, a = A[a], v = V[v],
                               w = W[w], eps = err_tol), # only "lower" resp
          times = times, unit = unit)
        # add the v, a, w values and function names to the dataframe
        mbm_res[start:stop, 1] <- rep(V[v]  , nf)
        mbm_res[start:stop, 2] <- rep(A[a]  , nf)
        mbm_res[start:stop, 3] <- rep(W[w]  , nf)
        mbm_res[start:stop, 4] <- fnames
        # add the microbenchmark results to the dataframe
        for (i in 1:nf) {
          bm <- subset(mbm, expr==fnames[i])
          mbm_res[start+i-1, sl_time] <- bm$time
        }
        # iterate start and stop values
        start = start + nf
        stop = stop + nf
      }
    }
  }
  return(mbm_res)
}










####################### Wrapper for Benchmarks #################################
rt_benchmark <- function(RT, resp, V, A, W, t0 = 1e-4, err_tol = 1e-6,
                         as_vec = FALSE, times = 1000, unit = "us")
{
  if (as_vec) { # benchmark all rt's as a vector
    return(bm_vec(RT = RT, resp = resp, V = V, A = A, W = W, t0 = t0,
                  err_tol = err_tol, times = times, unit = unit))
  } else { # benchmark each rt individually
    return(bm(RT = RT, resp = resp, V = V, A = A, W = W, t0 = t0,
              err_tol = err_tol, times = times, unit = unit))
  }
}
