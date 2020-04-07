
library("microbenchmark")
library("rtdists")
library("RWiener")

RT <- 1

A <- 1.5
a <- 1
V <- 2
v <- 1
t0 <- 1e-4 # must be nonzero for RWiener
W <- 0.5
w <- 1
err_tol = 1e-6 # this is the setting from rtdists
resp <- 0
rresp <- "lower"
source("inst/extdata/Kesselmeier_density.R")

mbm <- microbenchmark(
  fs_Fos_17 = dfddm(rt = RT, response = resp, a = A[a],
                    v = V[v], t0 = t0, w = W[w],
                    log = FALSE, n_terms_small = "Foster",
                    summation_small = "2017", scale = "small",
                    err_tol = err_tol),
  fb_Kes_17 = dfddm(rt = RT, response = resp, a = A[a],
                    v = V[v], t0 = t0, w = W[w],
                    log = FALSE, n_terms_small = "Kesselmeier",
                    summation_small = "2017", scale = "both",
                    err_tol = err_tol),
  fb_Nav_17 = dfddm(rt = RT, response = resp, a = A[a],
                    v = V[v], t0 = t0, w = W[w],
                    log = FALSE, n_terms_small = "Navarro",
                    summation_small = "2017", scale = "both",
                    err_tol = err_tol),
  RWiener = dwiener(RT, resp = rresp, alpha = A[a],
                    delta = V[v], tau = t0, beta = W[w],
                    give_log = FALSE),
  Kesselmeier = fs14_R(t = RT-t0, a = A[a], v = V[v],
                       w = W[w], eps = err_tol), check = "equal")
mbm
str(mbm)
