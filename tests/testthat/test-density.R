context("Testing the density function for accuracy")
### See Known Errors (KE) at bottom

library("rtdists")
library("RWiener")
source(system.file("extdata", "Gondan_et_al_density.R", package = "fddm", mustWork = TRUE))


### Evaluate densities for checking later ##
# Define different parameter spaces
if (identical(Sys.getenv("NOT_CRAN"), "true")) { # not on CRAN
  RT <- c(0.001, 0.1, 1, 2, 3, 4, 5, 10, 30)
  A <- c(0.25, 0.5, 1, 2.5, 5)
  V <- c(-5, -2, 0, 2, 5)
  W <- c(0.2, 0.5, 0.8)
  SV <- c(0, 0.5, 1, 1.5)
} else { # on CRAN
  RT <- c(0.001, 0.1, 1, 10)
  A <- c(0.5, 1, 5)
  V <- c(-5, 0, 5)
  W <- c(0.2, 0.5, 0.8)
  SV <- c(0, 0.5, 1.5)
}
t0 <- 1e-4 # must be nonzero for RWiener
SV_THRESH <- 1e-6
eps <- 1e-6 # this is the setting from rtdists

nRT <- length(RT)
nA <- length(A)
nV <- length(V)
nW <- length(W)
nSV <- length(SV)
resp <- rep("lower", nRT) # for RWiener

fnames <- c("fs_SWSE_17", "fs_SWSE_14", "fb_SWSE_17", "fb_SWSE_14",
            "fs_Gon_17", "fs_Gon_14", "fb_Gon_17", "fb_Gon_14",
            "fs_Nav_17", "fs_Nav_14", "fb_Nav_17", "fb_Nav_14",
            "fl_Nav_09", "RWiener", "Gondan", "rtdists")
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
          res[start,    7] <- dfddm(rt = RT[rt], response = resp[rt], a = A[a],
                                    v = V[v], t0 = t0, w = W[w], sv = SV[sv],
                                    log = FALSE, n_terms_small = "SWSE",
                                    summation_small = "2017", scale = "small",
                                    err_tol = eps)
          res[start+1,  7] <- dfddm(rt = RT[rt], response = resp[rt], a = A[a],
                                    v = V[v], t0 = t0, w = W[w], sv = SV[sv],
                                    log = FALSE, n_terms_small = "SWSE",
                                    summation_small = "2014", scale = "small",
                                    err_tol = eps)
          res[start+2,  7] <- dfddm(rt = RT[rt], response = resp[rt], a = A[a],
                                    v = V[v], t0 = t0, w = W[w], sv = SV[sv],
                                    log = FALSE, n_terms_small = "SWSE",
                                    summation_small = "2017", scale = "both",
                                    err_tol = eps)
          res[start+3,  7] <- dfddm(rt = RT[rt], response = resp[rt], a = A[a],
                                    v = V[v], t0 = t0, w = W[w], sv = SV[sv],
                                    log = FALSE, n_terms_small = "SWSE",
                                    summation_small = "2014", scale = "both",
                                    err_tol = eps)
          res[start+4,  7] <- dfddm(rt = RT[rt], response = resp[rt], a = A[a],
                                    v = V[v], t0 = t0, w = W[w], sv = SV[sv],
                                    log = FALSE, n_terms_small = "Gondan",
                                    summation_small = "2017", scale = "small",
                                    err_tol = eps)
          res[start+5,  7] <- dfddm(rt = RT[rt], response = resp[rt], a = A[a],
                                    v = V[v], t0 = t0, w = W[w], sv = SV[sv],
                                    log = FALSE, n_terms_small = "Gondan",
                                    summation_small = "2014", scale = "small",
                                    err_tol = eps)
          res[start+6,  7] <- dfddm(rt = RT[rt], response = resp[rt], a = A[a],
                                    v = V[v], t0 = t0, w = W[w], sv = SV[sv],
                                    log = FALSE, n_terms_small = "Gondan",
                                    summation_small = "2017", scale = "both",
                                    err_tol = eps)
          res[start+7,  7] <- dfddm(rt = RT[rt], response = resp[rt], a = A[a],
                                    v = V[v], t0 = t0, w = W[w], sv = SV[sv],
                                    log = FALSE, n_terms_small = "Gondan",
                                    summation_small = "2014", scale = "both",
                                    err_tol = eps)
          res[start+8,  7] <- dfddm(rt = RT[rt], response = resp[rt], a = A[a],
                                    v = V[v], t0 = t0, w = W[w], sv = SV[sv],
                                    log = FALSE, n_terms_small = "Navarro",
                                    summation_small = "2017", scale = "small",
                                    err_tol = eps)
          res[start+9,  7] <- dfddm(rt = RT[rt], response = resp[rt], a = A[a],
                                    v = V[v], t0 = t0, w = W[w], sv = SV[sv],
                                    log = FALSE, n_terms_small = "Navarro",
                                    summation_small = "2014", scale = "small",
                                    err_tol = eps)
          res[start+10, 7] <- dfddm(rt = RT[rt], response = resp[rt], a = A[a],
                                    v = V[v], t0 = t0, w = W[w], sv = SV[sv],
                                    log = FALSE, n_terms_small = "Navarro",
                                    summation_small = "2017", scale = "both",
                                    err_tol = eps)
          res[start+11, 7] <- dfddm(rt = RT[rt], response = resp[rt], a = A[a],
                                    v = V[v], t0 = t0, w = W[w], sv = SV[sv],
                                    log = FALSE, n_terms_small = "Navarro",
                                    summation_small = "2014", scale = "both",
                                    err_tol = eps)
          res[start+12,  7] <- dfddm(rt = RT[rt], response = resp[rt], a = A[a],
                                    v = V[v], t0 = t0, w = W[w], sv = SV[sv],
                                    log = FALSE, n_terms_small = "Navarro",
                                    scale = "large", err_tol = eps)
          res[start+13, 7] <- dwiener(RT[rt], resp = resp[rt], alpha = A[a],
                                      delta = V[v], tau = t0, beta = W[w],
                                      give_log = FALSE)
          res[start+14, 7] <- fs(t = RT[rt]-t0, a = A[a], v = V[v],
                                 w = W[w], eps = eps)
          res[start+15, 7] <- ddiffusion(RT[rt], resp[rt], a = A[a], v = V[v],
                                        t0 = t0, z = W[w]*A[a], sv = SV[sv])
          if (sv > SV_THRESH) { # multiply to get density with sv
            t <- RT[rt] - t0
            M <- exp(V[v] * A[a] * W[w] + V[v]*V[v] * t / 2 +
                     (SV[sv]*SV[sv] * A[a]*A[a] * W[w]*W[w] -
                       2 * V[v] * A[a] * W[w] - V[v]*V[v] * t) /
                     (2 + 2 * SV[sv]*SV[sv] * t)) / sqrt(1 + SV[sv]*SV[sv] * t)
            res[start+13, 7] <- M * res[start+11, 7] # RWiener
            res[start+14, 7] <- M * res[start+12, 7] # Gondan_R
          }

          # calculate differences
          ans <- res[start + 2, 7] # use fb_SWSE_17 as truth
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
          res[start+14, 8] <- abs(res[start+14, 7] - ans)
          res[start+15, 8] <- abs(res[start+15, 7] - ans)

          # calculate log of "lower" density
          res[start,    9] <- dfddm(rt = RT[rt], response = resp[rt], a = A[a],
                                    v = V[v], t0 = t0, w = W[w], sv = SV[sv],
                                    log = TRUE, n_terms_small = "SWSE",
                                    summation_small = "2017", scale = "small",
                                    err_tol = eps)
          res[start+1,  9] <- dfddm(rt = RT[rt], response = resp[rt], a = A[a],
                                    v = V[v], t0 = t0, w = W[w], sv = SV[sv],
                                    log = TRUE, n_terms_small = "SWSE",
                                    summation_small = "2014", scale = "small",
                                    err_tol = eps)
          res[start+2,  9] <- dfddm(rt = RT[rt], response = resp[rt], a = A[a],
                                    v = V[v], t0 = t0, w = W[w], sv = SV[sv],
                                    log = TRUE, n_terms_small = "SWSE",
                                    summation_small = "2017", scale = "both",
                                    err_tol = eps)
          res[start+3,  9] <- dfddm(rt = RT[rt], response = resp[rt], a = A[a],
                                    v = V[v], t0 = t0, w = W[w], sv = SV[sv],
                                    log = TRUE, n_terms_small = "SWSE",
                                    summation_small = "2014", scale = "both",
                                    err_tol = eps)
          res[start+4,  9] <- dfddm(rt = RT[rt], response = resp[rt], a = A[a],
                                    v = V[v], t0 = t0, w = W[w], sv = SV[sv],
                                    log = TRUE, n_terms_small = "Gondan",
                                    summation_small = "2017", scale = "small",
                                    err_tol = eps)
          res[start+5,  9] <- dfddm(rt = RT[rt], response = resp[rt], a = A[a],
                                    v = V[v], t0 = t0, w = W[w], sv = SV[sv],
                                    log = TRUE, n_terms_small = "Gondan",
                                    summation_small = "2014", scale = "small",
                                    err_tol = eps)
          res[start+6,  9] <- dfddm(rt = RT[rt], response = resp[rt], a = A[a],
                                    v = V[v], t0 = t0, w = W[w], sv = SV[sv],
                                    log = TRUE, n_terms_small = "Gondan",
                                    summation_small = "2017", scale = "both",
                                    err_tol = eps)
          res[start+7,  9] <- dfddm(rt = RT[rt], response = resp[rt], a = A[a],
                                    v = V[v], t0 = t0, w = W[w], sv = SV[sv],
                                    log = TRUE, n_terms_small = "Gondan",
                                    summation_small = "2014", scale = "both",
                                    err_tol = eps)
          res[start+8,  9] <- dfddm(rt = RT[rt], response = resp[rt], a = A[a],
                                    v = V[v], t0 = t0, w = W[w], sv = SV[sv],
                                    log = TRUE, n_terms_small = "Navarro",
                                    summation_small = "2017", scale = "small",
                                    err_tol = eps)
          res[start+9,  9] <- dfddm(rt = RT[rt], response = resp[rt], a = A[a],
                                    v = V[v], t0 = t0, w = W[w], sv = SV[sv],
                                    log = TRUE, n_terms_small = "Navarro",
                                    summation_small = "2014", scale = "small",
                                    err_tol = eps)
          res[start+10, 9] <- dfddm(rt = RT[rt], response = resp[rt], a = A[a],
                                    v = V[v], t0 = t0, w = W[w], sv = SV[sv],
                                    log = TRUE, n_terms_small = "Navarro",
                                    summation_small = "2017", scale = "both",
                                    err_tol = eps)
          res[start+11, 9] <- dfddm(rt = RT[rt], response = resp[rt], a = A[a],
                                    v = V[v], t0 = t0, w = W[w], sv = SV[sv],
                                    log = TRUE, n_terms_small = "Navarro",
                                    summation_small = "2014", scale = "both",
                                    err_tol = eps)
          res[start+12,  9] <- dfddm(rt = RT[rt], response = resp[rt], a = A[a],
                                    v = V[v], t0 = t0, w = W[w], sv = SV[sv],
                                    log = TRUE, n_terms_small = "Navarro",
                                    scale = "large", err_tol = eps)
          res[start+13, 9] <- dwiener(RT[rt], resp = resp[rt], alpha = A[a],
                                      delta = V[v], tau = t0, beta = W[w],
                                      give_log = TRUE)
          res[start+14, 9] <- log(fs(t = RT[rt]-t0, a = A[a], v = V[v],
                                     w = W[w], eps = eps))
          res[start+15, 9] <- log(ddiffusion(RT[rt], resp[rt], a = A[a],
                                             v = V[v], t0 = t0, z = W[w]*A[a],
                                             sv = SV[sv]))
          if (sv > SV_THRESH) { # add to get log of density with sv
            t <- RT[rt] - t0
            M <- V[v] * A[a] * W[w] + V[v]*V[v] * t / 2 +
                 (SV[sv]*SV[sv] * A[a]*A[a] * W[w]*W[w] -
                  2 * V[v] * A[a] * W[w] - V[v]*V[v] * t) /
                 (2 + 2 * SV[sv]*SV[sv] * t) - 0.5 * log(1 + SV[sv]*SV[sv] * t)
            res[start+13, 9] <- M + res[start+11, 9] # RWiener
            res[start+14, 9] <- M + res[start+12, 9] # Gondan_R
          }

          # iterate start and stop values
          start = start + nf
          stop = stop + nf
        }
      }
    }
  }
}


### Prep for testing ###
# Subset results
SWSE_s <- res[res[["FuncName"]] %in% fnames[c(1, 2)], ]
SWSE_b <- res[res[["FuncName"]] %in% fnames[c(3, 4)], ]
Gondan_s <- res[res[["FuncName"]] %in% fnames[c(5, 6)], ]
Gondan_b <- res[res[["FuncName"]] %in% fnames[c(7, 8)], ]
Navarro_s <- res[res[["FuncName"]] %in% fnames[c(9, 10)], ]
Navarro_b <- res[res[["FuncName"]] %in% fnames[c(11, 12)], ]
Navarro_l <- res[res[["FuncName"]] %in% fnames[13], ]
RWiener <- res[res[["FuncName"]] %in% fnames[14], ]
Gondan_R <- res[res[["FuncName"]] %in% fnames[15], ]
rtdists <- res[res[["FuncName"]] %in% fnames[16], ]


### Testing ###
# Ensure all densities are non-negative
test_that("Non-negativity of densities", {
  expect_true(all(SWSE_s[["res"]] >= 0))
  expect_true(all(SWSE_b[["res"]] >= 0))
  expect_true(all(Gondan_s[["res"]] >= 0))
  expect_true(all(Gondan_b[["res"]] >= 0))
  expect_true(all(Navarro_s[["res"]] >= 0))
  expect_true(all(Navarro_b[["res"]] >= 0))
  expect_true(all(Navarro_l[["res"]] >= 0))
  expect_true(all(RWiener[["res"]] >= 0))
  expect_true(all(Gondan_R[["res"]] >= 0))
  expect_true(all(rtdists[["res"]] >= 0))
})

# Test accuracy within 2*eps (allows for convergence from above and below)
test_that("Consistency among internal methods", {
  expect_true(all(SWSE_s[["dif"]] < 2*eps))
  expect_true(all(SWSE_b[["dif"]] < 2*eps))
  expect_true(all(Gondan_s[["dif"]] < 2*eps))
  expect_true(all(Gondan_b[["dif"]] < 2*eps))
  expect_true(all(Navarro_s[["dif"]] < 2*eps))
  expect_true(all(Navarro_b[["dif"]] < 2*eps))
  expect_true(all(Navarro_l[["dif"]] < 2*eps))
})

test_that("Accuracy relative to established packages", {
  expect_true(all(RWiener[RWiener[["sv"]] < SV_THRESH, "dif"] < 2*eps)) # see KE 1
  expect_true(all(Gondan_R[Gondan_R[["sv"]] < SV_THRESH, "dif"] < 2*eps)) # see KE 1
  expect_true(all(rtdists[["dif"]] < 2*eps))
})

# Test consistency in log vs non-log (see KE 2)
test_that("Log-Consistency among internal methods", {
  expect_equal(SWSE_s[SWSE_s[["res"]] > eps*eps, "log_res"],
               log(SWSE_s[SWSE_s[["res"]] > eps*eps, "res"]))
  expect_equal(SWSE_b[SWSE_b[["res"]] > eps*eps, "log_res"],
               log(SWSE_b[SWSE_b[["res"]] > eps*eps, "res"]))
  expect_equal(Gondan_s[Gondan_s[["res"]] > eps*eps, "log_res"],
               log(Gondan_s[Gondan_s[["res"]] > eps*eps, "res"]))
  expect_equal(Gondan_b[Gondan_b[["res"]] > eps*eps, "log_res"],
               log(Gondan_b[Gondan_b[["res"]] > eps*eps, "res"]))
  expect_equal(Navarro_s[Navarro_s[["res"]] > eps*eps, "log_res"],
               log(Navarro_s[Navarro_s[["res"]] > eps*eps, "res"]))
  expect_equal(Navarro_b[Navarro_b[["res"]] > eps*eps, "log_res"],
               log(Navarro_b[Navarro_b[["res"]] > eps*eps, "res"]))
  expect_equal(Navarro_l[Navarro_l[["res"]] > eps*eps, "log_res"],
               log(Navarro_l[Navarro_l[["res"]] > eps*eps, "res"]))
})

test_that("Log-Consistency of established packages", {
  testthat::skip_on_cran()
  expect_equal(RWiener[RWiener[["res"]] > eps*eps, "log_res"],
               log(RWiener[RWiener[["res"]] > eps*eps, "res"]))
  expect_equal(Gondan_R[Gondan_R[["res"]] > eps*eps, "log_res"],
               log(Gondan_R[Gondan_R[["res"]] > eps*eps, "res"]))
  expect_equal(rtdists[rtdists[["res"]] > eps*eps, "log_res"],
               log(rtdists[rtdists[["res"]] > eps*eps, "res"]))
})



### Known Errors (KE) ###
#
# 1) Both RWiener and Gondan_R divide the error tolerance by the multiplicative
#    term outside of the summation. Since the outside term is different when
#    $sv > 0$, the approximations use the incorrect error tolerance for
#    $sv > 0$. This affects the number of terms required in the summation to
#    achieve the desired precision, thus not actually achieving that desired
#    precision. This issue is fixed in our implementation of the Gondan method,
#    `n_terms_small = "Gondan"`, `scale = "small"`. For an example of this
#    discrepancy, see the code below:
#
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
# ks(t/(a*a), w, eps/sv0) # = 2; the summation will only calculate 2 terms
# ks(t/(a*a), w, eps/sv0_9) # = 5; but the summation actually needs 5 terms
#
# 2) When calculating the log of the density, it is better to use the built-in
#    log option. For very small densities, simply calculating the density can
#    cause rounding issues that result in a density of zero (thus the log of the
#    density becomes -Inf). Using the built-in log option avoids some of these
#    rounding issues by exploiting the algebraic properties of the logarithm.
#    Also note that sometimes the densities are just too small (i.e. extremely
#    negative) and the logarithm function returns a value of -Inf, so we discard
#    the samples whose density is very small (less than eps*eps = 1e-12).
