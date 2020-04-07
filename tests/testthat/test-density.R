context("Accuracy among internal methods and relative to established packages")
library("rtdists")
library("RWiener")
source("inst/extdata/Kesselmeier_density.R")

### See Known Errors (KE) at bottom for issues with Navarro's implementation



# Evaluate densities
RT = c(0.001, 0.005, 0.01, 0.05, 0.1,
       seq(0.5, 3, by = 0.5),
       seq(4, 10, by = 1),
       seq(12.5, 20, by = 2.5),
       seq(25, 30, by = 5))
A = seq(0.5, 5, by = 0.5)
V = seq(-6, 6, by = 0.5)
t0 = 1e-4 # must be nonzero for RWiener
W = seq(0.2, 0.8, by = 0.1)
SV = seq(0, 3, by = 0.3)
eps = 1e-6 # this is the setting from rtdists

nRT <- length(RT)
nA <- length(A)
nV <- length(V)
nW <- length(W)
nSV <- length(SV)
resp0 = rep(0, nRT)
rresp0 = rep("lower", nRT)

nf <- 14
fnames <- c("fs_Fos_17", "fs_Fos_14", "fs_Kes_17", "fs_Kes_14",
            "fs_Nav_17", "fs_Nav_14", "fl_Nav_09",
            "fb_Kes_17", "fb_Kes_14", "fb_Nav_17", "fb_Nav_14",
            "rtdists", "RWiener", "Kesselmeier")

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
                                    log = FALSE, n_terms_small = NULL,
                                    summation_small = NULL, scale = "large",
                                    err_tol = eps)
          res[start+7,  7] <- dfddm(rt = RT[rt], response = resp0[rt], a = A[a],
                                    v = V[v], t0 = t0, w = W[w], sv = SV[sv],
                                    log = FALSE, n_terms_small = "Kesselmeier",
                                    summation_small = "2017", scale = "both",
                                    err_tol = eps)
          res[start+8,  7] <- dfddm(rt = RT[rt], response = resp0[rt], a = A[a],
                                    v = V[v], t0 = t0, w = W[w], sv = SV[sv],
                                    log = FALSE, n_terms_small = "Kesselmeier",
                                    summation_small = "2014", scale = "both",
                                    err_tol = eps)
          res[start+9,  7] <- dfddm(rt = RT[rt], response = resp0[rt], a = A[a],
                                    v = V[v], t0 = t0, w = W[w], sv = SV[sv],
                                    log = FALSE, n_terms_small = "Navarro",
                                    summation_small = "2017", scale = "both",
                                    err_tol = eps)
          res[start+10, 7] <- dfddm(rt = RT[rt], response = resp0[rt], a = A[a],
                                    v = V[v], t0 = t0, w = W[w], sv = SV[sv],
                                    log = FALSE, n_terms_small = "Navarro",
                                    summation_small = "2014", scale = "both",
                                    err_tol = eps)
          res[start+11, 7] <- ddiffusion(RT[rt], rresp0[rt], a = A[a], v = V[v],
                                         t0 = t0, z = W[w]*A[a], sv = SV[sv])
          # no sv handling for the following two densities
          res[start+12, 7] <- dwiener(RT[rt], resp = rresp0[rt], alpha = A[a],
                                      delta = V[v], tau = t0, beta = W[w],
                                      give_log = FALSE)
          res[start+13, 7] <- fs14_R(t = RT[rt]-t0, a = A[a], v = V[v],
                                     w = W[w], eps = eps)
          if (sv > 0.05) { # multiply to get density with sv
            t <- RT[rt] - t0
            M <- exp(V[v] * A[a] * W[w] + V[v]*V[v] * t / 2 +
                     (SV[sv]*SV[sv] * A[a]*A[a] * W[w]*W[w] -
                       2 * V[v] * A[a] * W[w] - V[v]*V[v] * t) /
                     (2 + 2 * SV[sv]*SV[sv] * t)) / sqrt(1 + SV[sv]*SV[sv] * t)
            res[start+12, 7] <- M * res[start+12, 7]
            res[start+13, 7] <- M * res[start+13, 7]
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
                                    log = TRUE, n_terms_small = NULL,
                                    summation_small = NULL, scale = "large",
                                    err_tol = eps)
          res[start+7,  9] <- dfddm(rt = RT[rt], response = resp0[rt], a = A[a],
                                    v = V[v], t0 = t0, w = W[w], sv = SV[sv],
                                    log = TRUE, n_terms_small = "Kesselmeier",
                                    summation_small = "2017", scale = "both",
                                    err_tol = eps)
          res[start+8,  9] <- dfddm(rt = RT[rt], response = resp0[rt], a = A[a],
                                    v = V[v], t0 = t0, w = W[w], sv = SV[sv],
                                    log = TRUE, n_terms_small = "Kesselmeier",
                                    summation_small = "2014", scale = "both",
                                    err_tol = eps)
          res[start+9,  9] <- dfddm(rt = RT[rt], response = resp0[rt], a = A[a],
                                    v = V[v], t0 = t0, w = W[w], sv = SV[sv],
                                    log = TRUE, n_terms_small = "Navarro",
                                    summation_small = "2017", scale = "both",
                                    err_tol = eps)
          res[start+10, 9] <- dfddm(rt = RT[rt], response = resp0[rt], a = A[a],
                                    v = V[v], t0 = t0, w = W[w], sv = SV[sv],
                                    log = TRUE, n_terms_small = "Navarro",
                                    summation_small = "2014", scale = "both",
                                    err_tol = eps)
          res[start+11, 9] <- log(ddiffusion(RT[rt], rresp0[rt], a = A[a],
                                             v = V[v], t0 = t0, z = W[w]*A[a],
                                             sv = SV[sv]))
          # no sv handling for the following two densities
          res[start+12, 9] <- dwiener(RT[rt], resp = rresp0[rt], alpha = A[a],
                                      delta = V[v], tau = t0, beta = W[w],
                                      give_log = TRUE)
          res[start+13, 9] <- log(fs14_R(t = RT[rt]-t0, a = A[a], v = V[v],
                                         w = W[w], eps = eps))
          if (sv > 0.05) { # add to get log of density with sv
            t <- RT[rt] - t0
            M <- V[v] * A[a] * W[w] + V[v]*V[v] * t / 2 +
                 (SV[sv]*SV[sv] * A[a]*A[a] * W[w]*W[w] -
                  2 * V[v] * A[a] * W[w] - V[v]*V[v] * t) /
                 (2 + 2 * SV[sv]*SV[sv] * t) - log(sqrt(1 + SV[sv]*SV[sv] * t))
            res[start+12, 9] <- M + res[start+12, 9]
            res[start+13, 9] <- M + res[start+13, 9]
          }

          # iterate start and stop values
          start = start + nf
          stop = stop + nf
        }
      }
    }
  }
}


# Subset
Foster <- subset(res, FuncName %in% fnames[1:2])
Kesselmeier <- subset(res, FuncName %in% fnames[c(3,4)])
Navarro <- subset(res, FuncName %in% fnames[c(5,6,7)])
Both_time_scales <- subset(res, FuncName %in% fnames[c(8,9,10,11)])
rtdists <- subset(res, FuncName == "rtdists")
RWiener <- subset(res, FuncName == "RWiener")
Kesselmeier_R <- subset(res, FuncName == "Kesselmeier")


# Test accuracy within 2*eps (allows for convergence from above and below)
test_that("Accuracy among internal methods", {
  expect_true(all(Foster$dif < 2*eps))
  expect_true(all(Kesselmeier$dif < 2*eps))
  # expect_true(all(Navarro$dif < 2*eps)) # see KE 1 and KE 2
  # expect_true(all(Both_time_scales$dif < 2*eps)) # see KE 2
})

test_that("Accuracy relative to established packages", {
  expect_true(all(rtdists$dif < 2*eps))
  expect_true(all(RWiener$dif < 2*eps))
  # expect_true(all(Kesselmeier_R$dif < 2*eps)) # see KE 3
})

test_that("Log-Accuracy among internal methods", {
  expect_equal(Foster$log_res, log(Foster$res))
  expect_equal(Kesselmeier$log_res, log(Kesselmeier$res))
  # expect_equal(Navarro$log_res, log(Navarro$res)) # see KE 1 and KE 2
  # expect_equal(Both_time_scales$log_res, log(Both_time_scales$res)) # see KE 2
})





### Known Errors (KE) ###
# 1) Navarro small time approximation gives slightly negative results for some
#    parameter values. These negative results are still within 2*eps of the
#    accepted results.
# 2) Navarro large time approximation is unstable for t <= 0.5
# 3) Kesselmeier approximation divides the error tolerance by the multiplicative
#    term outside of the summation. Since the outside term is different when
#    sv > 0, the approximation uses the incorrect error tolerance for sv > 0.
#    This affects the number of terms required in the summation to achieve the
#    desired precision, thus not actually achieving that desired precision. For
#    an example of this discrepancy, see the code below:
# rt <- 1.5
# t <- rt - 1e-4
# a <- 0.5
# v <- 4.5
# w <- 0.5
# sv <- 0.9
# sv0 <- exp(-v*a*w - v*v*t/2) / (a*a)
# sv0_9 <- exp((-2*v*a*w - v*v*t + sv*sv*a*a*w*w)/(2 + 2*sv*sv*t)) /
#          (a*a*sqrt(1+sv*sv*t))
# ks14_R(t/(a*a), w, eps/sv0) # = 2; the summation will only calculate 2 terms
# ks14_R(t/(a*a), w, eps/sv0_9) # = 5; but the summation actually needs 5 terms
