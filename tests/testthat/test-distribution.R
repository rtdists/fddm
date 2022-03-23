context("Testing the distribution function for accuracy")

source(system.file("extdata", "Blurton_et_al_distribution.R",
                   package = "fddm", mustWork = TRUE))
# note that pdiffusion() from rtdists is not agreeing with the Blurton/fddm code



#---------------------- Evaluate densities ------------------------------------#
# Define parameter space
if (identical(Sys.getenv("NOT_CRAN"), "true")) { # not on CRAN
  RT <- c(0.001, 0.01, seq(0.1, 10, by = 0.1), seq(15, 30, by = 5))
  A <- c(0.25, seq(0.5, 5, by = 0.5))
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
N <- nRT * nA * nV * nW * nSV

rt <- rep(RT, each = nSV * nW * nV * nA, times = 1)
a  <- rep(A,  each = nSV * nW * nV, times = nRT)
v  <- rep(V,  each = nSV * nW, times = nRT * nA)
w  <- rep(W,  each = nSV, times = nRT * nA * nV)
sv <- rep(SV, each = 1, times = nRT * nA * nV * nW)




# fddm methods
Mills <- data.frame(
  res = pfddm(rt = rt, response = "lower", a = a, v = v, t0 = t0, w = w,
              sv = sv, err_tol = eps, log = FALSE, method = "Mills"),
  dif = numeric(N),
  log_res = pfddm(rt = rt, response = "lower", a = a, v = v, t0 = t0, w = w,
                  sv = sv, err_tol = eps, log = TRUE, method = "Mills")
)
suppressWarnings( # warnings can be produced for large a, sv
NCDF <- data.frame(
  res = pfddm(rt = rt, response = "lower", a = a, v = v, t0 = t0, w = w,
              sv = sv, err_tol = eps, log = FALSE, method = "NCDF"),
  dif = numeric(N),
  log_res = pfddm(rt = rt, response = "lower", a = a, v = v, t0 = t0, w = w,
                  sv = sv, err_tol = eps, log = TRUE, method = "NCDF")
)
)

# non-fddm methods
Blurton <- data.frame(
  res = G_0(t = rt - t0, a = a, nu = v, w = w, eta2 = sv*sv, sigma2 = 1,
            eps = eps),
  dif = numeric(N)
)
# if (require("rtdists")) {
#   rtdists <- data.frame(
#     res = pdiffusion(rt, "lower", a = a, v = v, t0 = t0, z = w*a, sv = sv,
#                      precision = 3),
#     dif = numeric(N)
#   )
# }



# Calculate differences (use fddm's Mills method as truth)
ans <- Mills[["res"]]
Mills[["dif"]] <- abs(Mills[["res"]] - ans)
NCDF[["dif"]] <- abs(NCDF[["res"]] - ans)
Blurton[["dif"]] <- abs(Blurton[["res"]] - ans)
# if (require("rtdists")) {
#   rtdists[["dif"]] <- abs(rtdists[["res"]] - ans)
# }





#---------------------- Testing -----------------------------------------------#
# Ensure -eps = 0 <= CDF <= 1 = 1 + eps
test_that("CDF in [0, 1]", {
  expect_true(all(Mills[["res"]] >= 0 & Mills[["res"]] <= 1))
  expect_true(all(NCDF[["res"]] >= 0 & NCDF[["res"]] <= 1))

  expect_true(all(Blurton[["res"]] >= -eps & Blurton[["res"]] <= 1 + eps))

  # if (require("rtdists")) {
  #   expect_true(all(rtdists[["res"]] >= -eps & rtdists[["res"]] <= 1 + eps))
  # }
})

# Test accuracy within 2*eps (allows for convergence from above and below)
test_that("Accuracy relative to established packages", {
  expect_true(all(Mills[["dif"]] <= 2*eps))
  expect_true(sum(NCDF[["dif"]] <= 2*eps) / nrow(NCDF) > 0.9999) # 3 disagree
  expect_true(all(NCDF[["dif"]] <= 5*eps)) # they're still quite close

  expect_true(all(Blurton[["dif"]] <= 2*eps))

  # if (require("rtdists")) {
  #   expect_true(all(rtdists[["dif"]] <= 2*eps))
  # }
})

# Test consistency in fddm log vs non-log
test_that("Log-Consistency", {
  expect_equal(Mills[Mills[["res"]] > eps*eps, "log_res"],
               log(Mills[Mills[["res"]] > eps*eps, "res"]))
  expect_equal(NCDF[NCDF[["res"]] > eps*eps, "log_res"],
               log(NCDF[NCDF[["res"]] > eps*eps, "res"]))
})
