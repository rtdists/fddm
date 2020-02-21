context("Accuracy among internal methods and relative to established packages")
library("rtdists")
library("RWiener")
source("inst/extdata/Kesselmeier_density.R")



# evaluate densities
RT = 0.1 # c(0.1, seq(0.5, 3, by = 0.5), seq(4, 10, by = 1))#,
# seq(12.5, 20, by = 2.5), seq(25, 30, b= 5))
# seq(1, 3, by = 0.5)
A = 0.5 # seq(0.5, 5, by = 0.5)
# seq(3, 6, by = 0.5)
V = -2.5 # # seq(-6, 6, by = 0.5)
# seq(-6, 6, by = 3.5)
t0 = 1e-4 # must be nonzero for RWiener
W = 0.5# seq(0.2, 0.8, by = 0.1)
# seq(0.2, 0.8, by = 0.3)
SV = # seq(0, 3, by = 0.5)
  # seq(0.5, 4.5, by = 1)
  0.1
eps = 1e-6 # this is the setting from rtdists

nRT <- length(RT) # 5
nA <- length(A) # 2
nV <- length(V) # 4
nW <- length(W) # 3
nSV <- length(SV) # 3
resp0 = rep(0, nRT)
rresp0 = rep("lower", nRT)

nf <- 14
fnames <- c("fs_Fos_17", "fs_Fos_14", "fs_Kes_17", "fs_Kes_14",
            "fs_Nav_17", "fs_Nav_14", "fl_Nav_09",
            "fb_Kes_17", "fb_Kes_14", "fb_Nav_17", "fb_Nav_14",
            "rtdists", "RWiener", "Kesselmeier")

res <- data.frame(matrix(ncol=8, nrow=nf*nRT*nA*nV*nW*nSV))
colnames(res) <- c('rt', 'a', 'v', 'w', 'sv', 'FuncName', 'res', 'dif')
start <- 1
stop <- nf

# Loop through each combination of parameters and record microbenchmark results
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
          res[start,    7] <- dfddm(rt = RT[rt], resp = resp0[rt], a = A[a],
                                    v = V[v], t0 = t0, w = W[w], sv = SV[sv],
                                    log = TRUE, n_terms_small = "Foster",
                                    summation_small = "2017", scale = "small",
                                    err_tol = eps)
          res[start+1,  7] <- dfddm(rt = RT[rt], resp = resp0[rt], a = A[a],
                                    v = V[v], t0 = t0, w = W[w], sv = SV[sv],
                                    log = TRUE, n_terms_small = "Foster",
                                    summation_small = "2014", scale = "small",
                                    err_tol = eps)
          res[start+2,  7] <- dfddm(rt = RT[rt], resp = resp0[rt], a = A[a],
                                    v = V[v], t0 = t0, w = W[w], sv = SV[sv],
                                    log = TRUE, n_terms_small = "Kesselmeier",
                                    summation_small = "2017", scale = "small",
                                    err_tol = eps)
          res[start+3,  7] <- dfddm(rt = RT[rt], resp = resp0[rt], a = A[a],
                                    v = V[v], t0 = t0, w = W[w], sv = SV[sv],
                                    log = TRUE, n_terms_small = "Kesselmeier",
                                    summation_small = "2014", scale = "small",
                                    err_tol = eps)
          res[start+4,  7] <- dfddm(rt = RT[rt], resp = resp0[rt], a = A[a],
                                    v = V[v], t0 = t0, w = W[w], sv = SV[sv],
                                    log = TRUE, n_terms_small = "Navarro",
                                    summation_small = "2017", scale = "small",
                                    err_tol = eps)
          res[start+5,  7] <- dfddm(rt = RT[rt], resp = resp0[rt], a = A[a],
                                    v = V[v], t0 = t0, w = W[w], sv = SV[sv],
                                    log = TRUE, n_terms_small = "Navarro",
                                    summation_small = "2014", scale = "small",
                                    err_tol = eps)
          res[start+6,  7] <- dfddm(rt = RT[rt], resp = resp0[rt], a = A[a],
                                    v = V[v], t0 = t0, w = W[w], sv = SV[sv],
                                    log = TRUE, n_terms_small = NULL,
                                    summation_small = NULL, scale = "large",
                                    err_tol = eps)
          res[start+7,  7] <- dfddm(rt = RT[rt], resp = resp0[rt], a = A[a],
                                    v = V[v], t0 = t0, w = W[w], sv = SV[sv],
                                    log = TRUE, n_terms_small = "Kesselmeier",
                                    summation_small = "2017", scale = "both",
                                    err_tol = eps)
          res[start+8,  7] <- dfddm(rt = RT[rt], resp = resp0[rt], a = A[a],
                                    v = V[v], t0 = t0, w = W[w], sv = SV[sv],
                                    log = TRUE, n_terms_small = "Kesselmeier",
                                    summation_small = "2014", scale = "both",
                                    err_tol = eps)
          res[start+9,  7] <- dfddm(rt = RT[rt], resp = resp0[rt], a = A[a],
                                    v = V[v], t0 = t0, w = W[w], sv = SV[sv],
                                    log = TRUE, n_terms_small = "Navarro",
                                    summation_small = "2017", scale = "both",
                                    err_tol = eps)
          res[start+10, 7] <- dfddm(rt = RT[rt], resp = resp0[rt], a = A[a],
                                    v = V[v], t0 = t0, w = W[w], sv = SV[sv],
                                    log = TRUE, n_terms_small = "Navarro",
                                    summation_small = "2014", scale = "both",
                                    err_tol = eps)
          res[start+11, 7] <- log(ddiffusion(RT[rt], rresp0[rt], a = A[a],
                                             v = V[v], t0 = t0, z = W[w]*A[a],
                                             sv = SV[sv]*SV[sv]))
          # no sv for following densities
          res[start+12, 7] <- dwiener(RT[rt], resp = rresp0[rt], alpha = A[a],
                                      delta = V[v], tau = t0, beta = W[w],
                                      give_log = TRUE)
          res[start+13, 7] <- log(fs14_R(t = RT[rt]-t0, a = A[a], v = V[v],
                                     w = W[w], eps = eps)) # "lower" boundary
          if (sv > 0) { # add to get log of density with sv
            t <- RT[rt] - t0
            M <- log(sqrt(2 * pi / (1 + SV[sv] * t)) *
                     exp(V[v] * A[a] * W[w] + V[v] * V[v] * t / 2 +
                         (SV[sv] * A[a] * A[a] * W[w] * W[w] -
                          2 * V[v] * A[a] * W[w] - V[v] * V[v] * t)
                         / (2 + 2 * SV[sv] * t)))
            # res[start+11, 7] <- M + res[start+11, 7] # rtdists
            res[start+12, 7] <- M + res[start+12, 7]
            res[start+13, 7] <- M + res[start+13, 7]
          }
          
          # calculate differences
          ans <- res[start+11, 7] # use rtdists as truth
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
          
          # iterate start and stop values
          start = start + nf
          stop = stop + nf
        }
      }
    }
  }
}


# test accuracy within 2*eps (allows for convergence from above and below)
test_that("Accuracy among internal methods", {
  small_time <- subset(res, FuncName %in% fnames[1:6])
  expect_true(all(small_time$dif < 2*eps))
  large_time <- subset(res, FuncName == fnames[7])
  expect_true(all(large_time$dif < 2*eps))
  both_times <- subset(res, FuncName %in% fnames[8:11])
  expect_true(all(both_times$dif < 2*eps))
})

test_that("Accuracy relative to established packages", {
  rtdists_dif <- subset(res, FuncName == "rtdists")
  expect_true(all(rtdists_dif$dif < 2*eps))
  RWiener_dif <- subset(res, FuncName == "RWiener")
  expect_true(all(RWiener_dif$dif < 2*eps))
  Kesselmeier_dif <- subset(res, FuncName == "Kesselmeier")
  expect_true(all(RWiener_dif$dif < 2*eps))
})
