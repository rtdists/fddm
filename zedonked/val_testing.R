devtools::load_all(recompile = TRUE)
library("fddm")
library("rtdists")





RT <- c(0.001, 0.1, 1, 10, 100)
A <- c(0.5, 1, 5)
V <- c(-5, 0, 5)
t0 <- 1e-4 # must be nonzero for RWiener
W <- c(0.2, 0.5, 0.8)
SV <- c(0, 0.5, 1.5)
SV_THRESH <- 1e-6
eps <- 1e-6 # this is the setting from rtdists

nRT <- length(RT)
nA <- length(A)
nV <- length(V)
nW <- length(W)
nSV <- length(SV)
resp <- rep("lower", nRT) # for RWiener

fnames <- c("fs_SWSE_17", "fb_Gon_17", "rtdists")
nf <- length(fnames)

res <- data.frame(matrix(ncol = 9, nrow = nf*nRT*nA*nV*nW*nSV))
colnames(res) <- c('rt', 'a', 'v', 'w', 'sv', 'FuncName', 'res', 'dif')
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
                                    log = FALSE, n_terms_small = "Gondan",
                                    summation_small = "2017", scale = "both",
                                    err_tol = eps)
          res[start+2, 7] <- ddiffusion(RT[rt], resp[rt], a = A[a], v = V[v],
                                        t0 = t0, z = W[w]*A[a], sv = SV[sv])

          # calculate differences
          ans <- res[start, 7] # use Foster_2017_small as truth
          res[start,    8] <- abs(res[start,    7] - ans)
          res[start+1,  8] <- abs(res[start+1,  7] - ans)
          res[start+2,  8] <- abs(res[start+2,  7] - ans)

          # iterate start and stop values
          start = start + nf
          stop = stop + nf
        }
      }
    }
  }
}


library("testthat")

SWSE <- subset(res, FuncName == "fc_SWSE_17")
Gondan_b <- subset(res, FuncName == "fb_Gon_17")
rtdists <- subset(res, FuncName == "rtdists")


expect_true(all(SWSE$res >= 0))
expect_true(all(Gondan_b$res >= 0))
expect_true(all(rtdists$res >= 0))

expect_true(all(SWSE$dif < 2*eps))
expect_true(all(Gondan_b$dif < 2*eps))
expect_true(all(rtdists$dif < 2*eps))

max(abs(Gondan_b$dif))
max(abs(rtdists$dif))
