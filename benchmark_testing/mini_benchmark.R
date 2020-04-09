
library("microbenchmark")
library("rtdists")
library("RWiener")
source("inst/extdata/Kesselmeier_density.R")
library("devtools")
load_all(recompile = TRUE)
library("fddm")

RT <- 1
A <- 1.5
V <- 2
t0 <- 1e-4 # must be nonzero for RWiener
W <- 0.5
err_tol = 1e-6 # this is the setting from rtdists
resp <- 0
rresp <- "lower"



dfddm(rt = RT, response = resp, a = A,
                  v = V, t0 = t0, w = W,
                  log = FALSE, n_terms_small = "Foster",
                  summation_small = "2017", scale = "small",
                  err_tol = err_tol)
dfddm_fast(rt = RT, a = A, v = V, t0 = t0, w = W, err_tol = err_tol)



mbm <- microbenchmark(
  fddm_fast = dfddm_fast(rt = RT, a = A, v = V, t0 = t0, w = W,
                         err_tol = err_tol),
  fs_Fos_17 = dfddm(rt = RT, response = resp, a = A,
                    v = V, t0 = t0, w = W,
                    log = FALSE, n_terms_small = "Foster",
                    summation_small = "2017", scale = "small",
                    err_tol = err_tol),
  fb_Kes_17 = dfddm(rt = RT, response = resp, a = A,
                    v = V, t0 = t0, w = W,
                    log = FALSE, n_terms_small = "Kesselmeier",
                    summation_small = "2017", scale = "both",
                    err_tol = err_tol),
  fb_Nav_17 = dfddm(rt = RT, response = resp, a = A,
                    v = V, t0 = t0, w = W,
                    log = FALSE, n_terms_small = "Navarro",
                    summation_small = "2017", scale = "both",
                    err_tol = err_tol),
  RWiener = dwiener(RT, resp = rresp, alpha = A,
                    delta = V, tau = t0, beta = W,
                    give_log = FALSE),
  Kesselmeier = fs14_R(t = RT-t0, a = A, v = V,
                       w = W, eps = err_tol),
  check = "equal")
summary(mbm)
