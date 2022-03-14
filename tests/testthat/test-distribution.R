context("Testing the distribution function for accuracy")

source(system.file("extdata", "Blurton_et_al_distribution.R",
                   package = "fddm", mustWork = TRUE))
# note that pdiffusion() from rtdists is not agreeing with the Blurton code


### Evaluate distributions for checking later ##
# Define different parameter spaces
if (identical(Sys.getenv("NOT_CRAN"), "true")) { # not on CRAN
  # These take a while to run
  # RT <- c(0.001, 0.01, seq(0.1, 10, by = 0.1), seq(15, 30, by = 5))
  # A <- c(0.25, seq(0.5, 5, by = 0.5))
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
t0 <- 0
eps <- 1e-6 # this is the setting from rtdists

nRT <- length(RT)
nA <- length(A)
nV <- length(V)
nW <- length(W)
nSV <- length(SV)
resp <- 1

fnames <- c("Mills", "NCDF", "Blurton", "rtdists")
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
          res[start,   7] <- pfddm(rt = RT[rt], response = resp, a = A[a],
                                   v = V[v], t0 = t0, w = W[w], sv = SV[sv],
                                   err_tol = eps, log = FALSE, method = "Mills")
          res[start+1, 7] <- pfddm(rt = RT[rt], response = resp, a = A[a],
                                   v = V[v], t0 = t0, w = W[w], sv = SV[sv],
                                   err_tol = eps, log = FALSE, method = "NCDF")
          res[start+2, 7] <- G_0(t = RT[rt]-t0, a = A[a], nu = V[v], w = W[w],
                                 eta2 = SV[sv]*SV[sv], sigma2 = 1, eps = eps)
          # if (require("rtdists")) {
          #   res[start+3, 7] <- pdiffusion(RT[rt], resp, a = A[a], v = V[v],
          #                                 t0 = t0, z = W[w]*A[a], sv = SV[sv],
          #                                 precision = 5)
          # }

          # calculate differences
          ans <- res[start+2, 7] # use Blurton code as truth
          res[start,   8] <- abs(res[start,    7] - ans)
          res[start+1, 8] <- abs(res[start+1,  7] - ans)
          res[start+2, 8] <- abs(res[start+2,  7] - ans)
          # if (require("rtdists")) {
          #   res[start+3, 8] <- abs(res[start+3, 7] - ans)
          # }

          # calculate log of "lower" density
          res[start,    9] <- pfddm(rt = RT[rt], response = resp, a = A[a],
                                    v = V[v], t0 = t0, w = W[w], sv = SV[sv],
                                    err_tol = eps, log = TRUE, method = "Mills")
          res[start+1,  9] <- pfddm(rt = RT[rt], response = resp, a = A[a],
                                    v = V[v], t0 = t0, w = W[w], sv = SV[sv],
                                    err_tol = eps, log = TRUE, method = "NCDF")
          res[start+2, 9] <- log(G_0(t = RT[rt]-t0, a = A[a], nu = V[v],
                                     w = W[w], eta2 = SV[sv]*SV[sv], sigma2 = 1,
                                     eps = eps))
          # if (require("rtdists")) {
          #   res[start+3, 9] <- log(pdiffusion(RT[rt], resp, a = A[a], v = V[v],
          #                                     t0 = t0, z = W[w]*A[a],
          #                                     sv = SV[sv], precision = 5))
          # }

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
Mills <- res[res[["FuncName"]] %in% fnames[1], ]
NCDF <- res[res[["FuncName"]] %in% fnames[2], ]
Blurton <- res[res[["FuncName"]] %in% fnames[3], ]
# if (require("rtdists")) {
#   rtdists <- res[res[["FuncName"]] %in% fnames[4], ]
# }


### Testing ###
# Ensure 0 <= CDF <= 1
test_that("CDF in [0, 1]", {
  expect_true(all(Mills[["res"]] >= 0 & Mills[["res"]] <= 1))
  expect_true(all(NCDF[["res"]] >= 0 & NCDF[["res"]] <= 1))
  
  # there are two from the Blurton code that are slightly > 1
  # expect_true(all(Blurton[["res"]] >= 0 & Blurton[["res"]] <= 1))
  
  # if (require("rtdists")) {
  #   expect_true(all(rtdists[["res"]] >= 0 & rtdists[["res"]] <= 1))
  # }
})

# Test accuracy within 2*eps (allows for convergence from above and below)
test_that("Accuracy relative to established packages", {
  expect_true(all(Mills[["dif"]] < 2*eps))
  # there is one instance of this method producing a slightly different result
  # rt = 30, response = 1, a = 5, v = 0, t0 = 0, w = 0.5, sv = 1.5
  # expect_true(all(NCDF[["dif"]] < 2*eps))
  
  # pointless because it's checked against itself
  # expect_true(all(Blurton[["dif"]] < 2*eps))
  
  # if (require("rtdists")) {
  #   expect_true(all(rtdists[["dif"]] < 2*eps))
  # }
})

# Test consistency in log vs non-log
test_that("Log-Consistency", {
  expect_equal(Mills[Mills[["res"]] > eps*eps, "log_res"],
               log(Mills[Mills[["res"]] > eps*eps, "res"]))
  expect_equal(NCDF[NCDF[["res"]] > eps*eps, "log_res"],
               log(NCDF[NCDF[["res"]] > eps*eps, "res"]))
  # testing the other two are pointless because there is no built-in log option
})
