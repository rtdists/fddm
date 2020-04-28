context("Accuracy among internal methods and relative to established packages")
library("rtdists")
library("RWiener")
source("inst/extdata/Kesselmeier_density.R")
library("devtools")
load_all(recompile = TRUE)
library("fddm")

### See Known Errors (KE) at bottom for issues with Navarro's implementation
### Can read result dataframe in following line or rerun in lines 14-232
# res <- readRDS("inst/extdata/test-density-res.Rds")


# Evaluate densities
RT <- c(0.001, 0.005, 0.01, 0.05, 0.1,
        seq(0.5, 3, by = 0.5),
        seq(4, 10, by = 1),
        seq(12.5, 20, by = 2.5),
        seq(25, 30, by = 5))
A <- seq(0.5, 5, by = 0.5)
V <- seq(-6, 6, by = 0.5)
t0 <- 1e-4 # must be nonzero for RWiener
W <- seq(0.2, 0.8, by = 0.1)
SV <- seq(0, 3, by = 0.3)
SV_THRESH <- 0.05
eps <- 1e-6 # this is the setting from rtdists

nRT <- length(RT)
nA <- length(A)
nV <- length(V)
nW <- length(W)
nSV <- length(SV)
resp0 <- rep(0, nRT)
rresp0 <- rep("lower", nRT)

fnames <- c("fs_Fos_17", "fs_Fos_14",
            "fs_Kes_17", "fs_Kes_14", "fs_Nav_17", "fs_Nav_14",
            "fb_Kes_17", "fb_Kes_14", "fb_Nav_17", "fb_Nav_14",
            "fl_Nav_09", "RWiener", "Kesselmeier", "rtdists")
nf <- length(fnames)

res <- data.frame(matrix(ncol = 9, nrow = nf*nRT*nA*nV*nW*nSV))
colnames(res) <- c('rt', 'a', 'v', 'w', 'sv', 'FuncName', 'res', 'dif',
                   'log_res')
start <- 1
stop <- nf

# Loop through each combination of parameters and record results
for (rt in 1:nRT) {
  for (a in 1:nA) {
    for (v in 1:nV) {
      for (w in 1:nW) {
        for (sv in 1:nSV) {
          # add the rt, v, a, w, and function names to the dataframe
          res[start:stop, 1] <- rep(RT[rt], nf)
          res[start:stop, 2] <- rep(A[a]  , nf)
          res[start:stop, 3] <- rep(V[v]  , nf)
          res[start:stop, 4] <- rep(W[w]  , nf)
          res[start:stop, 5] <- rep(SV[sv], nf)
          res[start:stop, 6] <- fnames

          # calculate "lower" density
          res[start,    7] <- dfddm(rt = RT[rt], response = resp0[rt], a = A[a],
                                    v = V[v], t0 = t0, w = W[w], sv = SV[sv],
                                    log = FALSE, n_terms_small = "Foster",
                                    summation_small = "2017", scale = "small",
                                    err_tol = eps)
          res[start+1,  7] <- dfddm(rt = RT[rt], response = resp0[rt], a = A[a],
                                    v = V[v], t0 = t0, w = W[w], sv = SV[sv],
                                    log = FALSE, n_terms_small = "Foster",
                                    summation_small = "2014", scale = "small",
                                    err_tol = eps)
          res[start+2,  7] <- dfddm(rt = RT[rt], response = resp0[rt], a = A[a],
                                    v = V[v], t0 = t0, w = W[w], sv = SV[sv],
                                    log = FALSE, n_terms_small = "Kesselmeier",
                                    summation_small = "2017", scale = "small",
                                    err_tol = eps)
          res[start+3,  7] <- dfddm(rt = RT[rt], response = resp0[rt], a = A[a],
                                    v = V[v], t0 = t0, w = W[w], sv = SV[sv],
                                    log = FALSE, n_terms_small = "Kesselmeier",
                                    summation_small = "2014", scale = "small",
                                    err_tol = eps)
          res[start+4,  7] <- dfddm(rt = RT[rt], response = resp0[rt], a = A[a],
                                    v = V[v], t0 = t0, w = W[w], sv = SV[sv],
                                    log = FALSE, n_terms_small = "Navarro",
                                    summation_small = "2017", scale = "small",
                                    err_tol = eps)
          res[start+5,  7] <- dfddm(rt = RT[rt], response = resp0[rt], a = A[a],
                                    v = V[v], t0 = t0, w = W[w], sv = SV[sv],
                                    log = FALSE, n_terms_small = "Navarro",
                                    summation_small = "2014", scale = "small",
                                    err_tol = eps)
          res[start+6,  7] <- dfddm(rt = RT[rt], response = resp0[rt], a = A[a],
                                    v = V[v], t0 = t0, w = W[w], sv = SV[sv],
                                    log = FALSE, n_terms_small = "Kesselmeier",
                                    summation_small = "2017", scale = "both",
                                    err_tol = eps)
          res[start+7,  7] <- dfddm(rt = RT[rt], response = resp0[rt], a = A[a],
                                    v = V[v], t0 = t0, w = W[w], sv = SV[sv],
                                    log = FALSE, n_terms_small = "Kesselmeier",
                                    summation_small = "2014", scale = "both",
                                    err_tol = eps)
          res[start+8,  7] <- dfddm(rt = RT[rt], response = resp0[rt], a = A[a],
                                    v = V[v], t0 = t0, w = W[w], sv = SV[sv],
                                    log = FALSE, n_terms_small = "Navarro",
                                    summation_small = "2017", scale = "both",
                                    err_tol = eps)
          res[start+9,  7] <- dfddm(rt = RT[rt], response = resp0[rt], a = A[a],
                                    v = V[v], t0 = t0, w = W[w], sv = SV[sv],
                                    log = FALSE, n_terms_small = "Navarro",
                                    summation_small = "2014", scale = "both",
                                    err_tol = eps)
          res[start+10,  7] <- dfddm(rt = RT[rt], response = resp0[rt], a = A[a],
                                    v = V[v], t0 = t0, w = W[w], sv = SV[sv],
                                    log = FALSE, n_terms_small = "",
                                    summation_small = "", scale = "large",
                                    err_tol = eps)
          res[start+11, 7] <- dwiener(RT[rt], resp = rresp0[rt], alpha = A[a],
                                      delta = V[v], tau = t0, beta = W[w],
                                      give_log = FALSE)
          res[start+12, 7] <- fs14_R(t = RT[rt]-t0, a = A[a], v = V[v],
                                     w = W[w], eps = eps)
          res[start+13, 7] <- ddiffusion(RT[rt], rresp0[rt], a = A[a], v = V[v],
                                        t0 = t0, z = W[w]*A[a], sv = SV[sv])
          if (sv > SV_THRESH) { # multiply to get density with sv
            t <- RT[rt] - t0
            M <- exp(V[v] * A[a] * W[w] + V[v]*V[v] * t / 2 +
                     (SV[sv]*SV[sv] * A[a]*A[a] * W[w]*W[w] -
                       2 * V[v] * A[a] * W[w] - V[v]*V[v] * t) /
                     (2 + 2 * SV[sv]*SV[sv] * t)) / sqrt(1 + SV[sv]*SV[sv] * t)
            res[start+11, 7] <- M * res[start+11, 7] # RWiener
            res[start+12, 7] <- M * res[start+12, 7] # Kesselmeier_R
          }

          # calculate differences
          ans <- res[start, 7] # use Foster_2017_small as truth
          res[start,    8] <- abs(res[start,    7] - ans)
          res[start+1,  8] <- abs(res[start+1,  7] - ans)
          res[start+2,  8] <- abs(res[start+2,  7] - ans)
          res[start+3,  8] <- abs(res[start+3,  7] - ans)
          res[start+4,  8] <- abs(res[start+4,  7] - ans)
          res[start+5,  8] <- abs(res[start+1,  7] - ans)
          res[start+6,  8] <- abs(res[start+6,  7] - ans)
          res[start+7,  8] <- abs(res[start+7,  7] - ans)
          res[start+8,  8] <- abs(res[start+8,  7] - ans)
          res[start+9,  8] <- abs(res[start+9,  7] - ans)
          res[start+10, 8] <- abs(res[start+10, 7] - ans)
          res[start+11, 8] <- abs(res[start+11, 7] - ans)
          res[start+12, 8] <- abs(res[start+12, 7] - ans)
          res[start+13, 8] <- abs(res[start+13, 7] - ans)

          # calculate log of "lower" density
          res[start,    9] <- dfddm(rt = RT[rt], response = resp0[rt], a = A[a],
                                    v = V[v], t0 = t0, w = W[w], sv = SV[sv],
                                    log = TRUE, n_terms_small = "Foster",
                                    summation_small = "2017", scale = "small",
                                    err_tol = eps)
          res[start+1,  9] <- dfddm(rt = RT[rt], response = resp0[rt], a = A[a],
                                    v = V[v], t0 = t0, w = W[w], sv = SV[sv],
                                    log = TRUE, n_terms_small = "Foster",
                                    summation_small = "2014", scale = "small",
                                    err_tol = eps)
          res[start+2,  9] <- dfddm(rt = RT[rt], response = resp0[rt], a = A[a],
                                    v = V[v], t0 = t0, w = W[w], sv = SV[sv],
                                    log = TRUE, n_terms_small = "Kesselmeier",
                                    summation_small = "2017", scale = "small",
                                    err_tol = eps)
          res[start+3,  9] <- dfddm(rt = RT[rt], response = resp0[rt], a = A[a],
                                    v = V[v], t0 = t0, w = W[w], sv = SV[sv],
                                    log = TRUE, n_terms_small = "Kesselmeier",
                                    summation_small = "2014", scale = "small",
                                    err_tol = eps)
          res[start+4,  9] <- dfddm(rt = RT[rt], response = resp0[rt], a = A[a],
                                    v = V[v], t0 = t0, w = W[w], sv = SV[sv],
                                    log = TRUE, n_terms_small = "Navarro",
                                    summation_small = "2017", scale = "small",
                                    err_tol = eps)
          res[start+5,  9] <- dfddm(rt = RT[rt], response = resp0[rt], a = A[a],
                                    v = V[v], t0 = t0, w = W[w], sv = SV[sv],
                                    log = TRUE, n_terms_small = "Navarro",
                                    summation_small = "2014", scale = "small",
                                    err_tol = eps)
          res[start+6,  9] <- dfddm(rt = RT[rt], response = resp0[rt], a = A[a],
                                    v = V[v], t0 = t0, w = W[w], sv = SV[sv],
                                    log = TRUE, n_terms_small = "Kesselmeier",
                                    summation_small = "2017", scale = "both",
                                    err_tol = eps)
          res[start+7,  9] <- dfddm(rt = RT[rt], response = resp0[rt], a = A[a],
                                    v = V[v], t0 = t0, w = W[w], sv = SV[sv],
                                    log = TRUE, n_terms_small = "Kesselmeier",
                                    summation_small = "2014", scale = "both",
                                    err_tol = eps)
          res[start+8,  9] <- dfddm(rt = RT[rt], response = resp0[rt], a = A[a],
                                    v = V[v], t0 = t0, w = W[w], sv = SV[sv],
                                    log = TRUE, n_terms_small = "Navarro",
                                    summation_small = "2017", scale = "both",
                                    err_tol = eps)
          res[start+9,  9] <- dfddm(rt = RT[rt], response = resp0[rt], a = A[a],
                                    v = V[v], t0 = t0, w = W[w], sv = SV[sv],
                                    log = TRUE, n_terms_small = "Navarro",
                                    summation_small = "2014", scale = "both",
                                    err_tol = eps)
          res[start+10,  9] <- dfddm(rt = RT[rt], response = resp0[rt], a = A[a],
                                    v = V[v], t0 = t0, w = W[w], sv = SV[sv],
                                    log = TRUE, n_terms_small = "",
                                    summation_small = "", scale = "large",
                                    err_tol = eps)
          res[start+11, 9] <- dwiener(RT[rt], resp = rresp0[rt], alpha = A[a],
                                      delta = V[v], tau = t0, beta = W[w],
                                      give_log = TRUE)
          res[start+12, 9] <- log(fs14_R(t = RT[rt]-t0, a = A[a], v = V[v],
                                         w = W[w], eps = eps))
          res[start+13, 9] <- log(ddiffusion(RT[rt], rresp0[rt], a = A[a],
                                             v = V[v], t0 = t0, z = W[w]*A[a],
                                             sv = SV[sv]))
          if (sv > SV_THRESH) { # add to get log of density with sv
            t <- RT[rt] - t0
            M <- V[v] * A[a] * W[w] + V[v]*V[v] * t / 2 +
                 (SV[sv]*SV[sv] * A[a]*A[a] * W[w]*W[w] -
                  2 * V[v] * A[a] * W[w] - V[v]*V[v] * t) /
                 (2 + 2 * SV[sv]*SV[sv] * t) - 0.5 * log(1 + SV[sv]*SV[sv] * t)
            res[start+11, 9] <- M + res[start+11, 9] # RWiener
            res[start+12, 9] <- M + res[start+12, 9] # Kesselmeier_R
          }

          # iterate start and stop values
          start = start + nf
          stop = stop + nf
        }
      }
    }
  }
}
# saveRDS(res, file = "inst/extdata/test-density-res.Rds")


# Subset results
Foster <- subset(res, FuncName %in% fnames[c(1,2)])
Kesselmeier_s <- subset(res, FuncName %in% fnames[c(3,4)])
Kesselmeier_b <- subset(res, FuncName %in% fnames[c(7,8)])
Navarro_s <- subset(res, FuncName %in% fnames[c(5,6)])
Navarro_l <- subset(res, FuncName %in% fnames[11])
Navarro_b <- subset(res, FuncName %in% fnames[c(9,10)])
rtdists <- subset(res, FuncName == "rtdists")
RWiener <- subset(res, FuncName == "RWiener")
Kesselmeier_R <- subset(res, FuncName == "Kesselmeier")

# Compensate for KE 1, 2
Nav_s_res_0 <- subset(Navarro_s, res < 0) # KE 2
Nav_s_0 <- min(Nav_s_res_0$rt / Nav_s_res_0$a / Nav_s_res_0$a) # = 10
Nav_l_res_0 <- subset(Navarro_l, res < 0) # KE 3
Nav_l_0 <- max(Nav_l_res_0$rt / Nav_l_res_0$a / Nav_l_res_0$a) # = 0.06
Nav_l_dif_2eps <- subset(Navarro_l, dif > 2*eps) # KE 3
Nav_l_2 <- max(Nav_l_dif_2eps$rt / Nav_l_dif_2eps$a / Nav_l_dif_2eps$a) # = 0.008

# Ensure all densities are non-negative
test_that("Non-negativity of densities", {
  expect_true(all(Foster$res >= 0))
  expect_true(all(Kesselmeier_s$res >= 0))
  expect_true(all(Kesselmeier_b$res >= 0))
  expect_true(all(subset(Navarro_s, rt/a/a < Nav_s_0)$res >= 0)) # see KE 1
  expect_true(all(subset(Navarro_l, rt/a/a > Nav_l_0)$res >= 0)) # see KE 2
  expect_true(all(subset(Navarro_b, rt/a/a > Nav_l_0)$res >= 0)) # see KE 3
  expect_true(all(rtdists$res >= 0))
  expect_true(all(RWiener$res >= 0))
  expect_true(all(Kesselmeier_R$res >= 0))
})

# Test accuracy within 2*eps (allows for convergence from above and below)
test_that("Accuracy among internal methods", {
  expect_true(all(Foster$dif < 2*eps))
  expect_true(all(Kesselmeier_s$dif < 2*eps))
  expect_true(all(Kesselmeier_b$dif < 2*eps))
  expect_true(all(Navarro_s$dif < 2*eps)) # see KE 1
  expect_true(all(subset(Navarro_l, rt/a/a > Nav_l_2)$dif < 2*eps)) # see KE 2
  expect_true(all(Navarro_b$dif < 2*eps)) # see KE 1,2
})

test_that("Accuracy relative to established packages", {
  expect_true(all(rtdists$dif < 2*eps))
  expect_true(all(RWiener$dif < 2*eps))
  expect_true(all(subset(Kesselmeier_R, sv < SV_THRESH)$dif < 2*eps)) # see KE 4
})

test_that("Log-Consistency among internal methods", {
  expect_equal(Foster$log_res, log(Foster$res))
  expect_equal(Kesselmeier_s$log_res, log(Kesselmeier_s$res))
  expect_equal(Kesselmeier_b$log_res, log(Kesselmeier_b$res))
  expect_equal(subset(Navarro_s, rt/a/a < Nav_s_0)$log_res,
               log(subset(Navarro_s, rt/a/a < Nav_s_0)$res)) # see KE 1
  expect_equal(subset(Navarro_l, rt/a/a > Nav_l_0)$log_res,
               log(subset(Navarro_l, rt/a/a > Nav_l_0)$res)) # see KE 2
  expect_equal(subset(Navarro_b, rt/a/a > Nav_l_0)$log_res,
               log(subset(Navarro_b, rt/a/a > Nav_l_0)$res)) # see KE 1,2
})

test_that("Log-Consistency of established packages", {
  expect_equal(rtdists$log_res, log(rtdists$res))
  expect_equal(RWiener$log_res, log(RWiener$res))
  expect_equal(Kesselmeier_R$log_res, log(Kesselmeier_R$res))
})





### Known Errors (KE) ###
# 1) Navarro small time approximation is unstable for "large" effective response
#    times, rt/(*a*a) >= 10. It gives slightly negative densities for effective
#    response times in this range. These parameter values should not have an
#    effect on the "both" time scale because the large time approximation should
#    handle such locations in the parameter space. The negative results,
#    however, are still within 2*eps of the accepted results (basically zero).
# 2) Navarro large time approximation is unstable for "small" effective response
#    times, rt/(a*a) <= 0.06. It gives slightly negative densities for effective
#    response times in this range and gives inaccurate densities for effective
#    response times of rt/(a*a) <= 0.008. These parameter values should not have
#    an effect on the "both" time scale because the small time approximation
#    should handle such locations in the parameter space.
# 3) The Navarro "both" time scale switches between the small time and large
#    time approximations. This method relies too much on the unstable large time
#    approximation, leading to its instability; thus we only need to subset to
#    correct for the large time.
# 4) Kesselmeier_R approximation divides the error tolerance by the
#    multiplicative term outside of the summation. Since the outside term is
#    different when sv > 0, the approximation uses the incorrect error tolerance
#    for sv > 0. This affects the number of terms required in the summation to
#    achieve the desired precision, thus not actually achieving that desired
#    precision. This is issue is fixed in our implementation of the Kesselmeier
#    method. For an example of this discrepancy, see the code below:
# rt <- 1.5
# t <- rt - 1e-4
# a <- 0.5
# v <- 4.5
# w <- 0.5
# eps <- 1e-6
# sv <- 0.9
# sv0 <- exp(-v*a*w - v*v*t/2) / (a*a) # for constant drift rate
# sv0_9 <- exp((-2*v*a*w - v*v*t + sv*sv*a*a*w*w)/(2 + 2*sv*sv*t)) /
#          (a*a*sqrt(1+sv*sv*t)) # for variable drift rate
# ks14_R(t/(a*a), w, eps/sv0) # = 2; the summation will only calculate 2 terms
# ks14_R(t/(a*a), w, eps/sv0_9) # = 5; but the summation actually needs 5 terms
