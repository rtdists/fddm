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

  fnames <- c("fddm_fast", "fs_Fos_17", "fs_Fos_14",
              "fs_Kes_17", "fs_Kes_14", "fs_Nav_17", "fs_Nav_14",
              "fb_Kes_17", "fb_Kes_14", "fb_Nav_17", "fb_Nav_14",
              "fl_Nav_09", "RWiener", "Kesselmeier", "rtdists")
  nf <- length(fnames) # number of functions being benchmarked
  nRT <- length(RT) # number of response times
  nV <- length(V) # number of drift rates
  nA <- length(A) # number of boundary separations
  nW <- length(W) # number of starting points

  # Initialize the dataframe to contain the microbenchmark results
  mbm_res <- data.frame(matrix(ncol = 4+nf, nrow=nRT*nV*nA*nW))
  colnames(mbm_res) <- c('RT', 'V', 'A', 'W', fnames)
  row_idx <- 1

  # Loop through each combination of parameters and record microbenchmark results
  for (rt in 1:nRT) {
    for (v in 1:nV) {
      for (a in 1:nA) {
        for (w in 1:nW) {
          mbm <- microbenchmark(
            fddm_fast = dfddm_fast(rt = RT[rt], a = A[a], v = V[v], t0 = t0,
                                   w = W[w], err_tol = err_tol),
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
            fl_Nav_09 = dfddm(rt = RT[rt], response = resp, a = A[a],
                              v = V[v], t0 = t0, w = W[w],
                              log = FALSE, n_terms_small = "",
                              summation_small = "", scale = "large",
                              err_tol = err_tol),
            RWiener = dwiener(RT[rt], resp = rresp, alpha = A[a],
                              delta = V[v], tau = t0, beta = W[w],
                              give_log = FALSE),
            Kesselmeier = fs14_R(t = RT[rt]-t0, a = A[a], v = V[v],
                                 w = W[w], eps = err_tol), # only "lower" resp
            rtdists = ddiffusion(RT[rt], rresp, a = A[a], v = V[v],
                                 t0 = t0, z = W[w]*A[a]),
            times = times, unit = unit)
          # add the rt, v, a, and w values to the dataframe
          mbm_res[row_idx, 1] <- RT[rt]
          mbm_res[row_idx, 2] <- V[v]
          mbm_res[row_idx, 3] <- A[a]
          mbm_res[row_idx, 4] <- W[w]
          # add the median microbenchmark results to the dataframe
          for (i in 1:nf) {
            mbm_res[row_idx, 4+i] <- median(subset(mbm, expr==fnames[i])$time)
          }
          # iterate start value
          row_idx = row_idx + 1
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

  fnames <- c("fddm_fast", "fs_Fos_17", "fs_Fos_14",
              "fs_Kes_17", "fs_Kes_14", "fs_Nav_17", "fs_Nav_14",
              "fb_Kes_17", "fb_Kes_14", "fb_Nav_17", "fb_Nav_14",
              "fl_Nav_09", "RWiener", "Kesselmeier", "rtdists")
  nf <- length(fnames) # number of functions being benchmarked
  nV <- length(V) # number of drift rates
  nA <- length(A) # number of boundary separations
  nW <- length(W) # number of starting points
  rresp <- rep(rresp, length(RT)) # for RWiener

  # Initialize the dataframe to contain the microbenchmark results
  mbm_res <- data.frame(matrix(ncol = 3+nf, nrow = nV*nA*nW))
  colnames(mbm_res) <- c('V', 'A', 'W', fnames)
  row_idx <- 1

  # Loop through each combination of parameters and record microbenchmark results
  for (v in 1:nV) {
    for (a in 1:nA) {
      for (w in 1:nW) {
        mbm <- microbenchmark(
          fddm_fast = dfddm_fast(rt = RT, a = A[a], v = V[v], t0 = t0,
                                 w = W[w], err_tol = err_tol),
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
          fl_Nav_09 = dfddm(rt = RT, response = resp, a = A[a],
                            v = V[v], t0 = t0, w = W[w],
                            log = FALSE, n_terms_small = "",
                            summation_small = "", scale = "large",
                            err_tol = err_tol),
          RWiener = dwiener(RT, resp = rresp, alpha = A[a],
                            delta = V[v], tau = t0, beta = W[w],
                            give_log = FALSE),
          Kesselmeier = fs14_R(t = RT-t0, a = A[a], v = V[v],
                               w = W[w], eps = err_tol), # only "lower" resp
          rtdists = ddiffusion(RT, rresp, a = A[a], v = V[v],
                               t0 = t0, z = W[w]*A[a]),
          times = times, unit = unit)
        # add the v, a, and w values to the dataframe
        mbm_res[row_idx, 1] <- V[v]
        mbm_res[row_idx, 2] <- A[a]
        mbm_res[row_idx, 3] <- W[w]
        # add the median microbenchmark results to the dataframe
        for (i in 1:nf) {
          mbm_res[row_idx, 3+i] <- median(subset(mbm, expr==fnames[i])$time)
        }
        # iterate start value
        row_idx = row_idx + 1
      }
    }
  }
  return(mbm_res)
}










####################### Wrapper for Benchmarks #################################
rt_benchmark <- function(RT, resp, V, A, W, t0 = 1e-4, err_tol = 1e-6,
                         rt_as_vec = FALSE, times = 1000, unit = "us")
{
  if (rt_as_vec) { # benchmark all rt's as a vector
    return(bm_vec(RT = RT, resp = resp, V = V, A = A, W = W, t0 = t0,
                  err_tol = err_tol, times = times, unit = unit))
  } else { # benchmark each rt individually
    return(bm(RT = RT, resp = resp, V = V, A = A, W = W, t0 = t0,
              err_tol = err_tol, times = times, unit = unit))
  }
}
