context("Testing the density function for accuracy")
### See Known Errors (KE) at bottom

source(system.file("extdata", "Gondan_et_al_density.R",
                   package = "fddm", mustWork = TRUE))



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
t0 <- 1e-4 # must be nonzero for RWiener
SV_THRESH <- 1e-6
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
SWSE_s_17 <- data.frame(
  res = dfddm(rt = rt, response = "lower", a = a, v = v, t0 = t0, w = w,
              sv = sv, err_tol = eps, log = FALSE, switch_mech = "small",
              n_terms_small = "SWSE", summation_small = "2017"),
  dif = numeric(N),
  log_res = dfddm(rt = rt, response = "lower", a = a, v = v, t0 = t0, w = w,
                  sv = sv, err_tol = eps, log = TRUE, switch_mech = "small",
                  n_terms_small = "SWSE", summation_small = "2017")
)
SWSE_s_14 <- data.frame(
  res = dfddm(rt = rt, response = "lower", a = a, v = v, t0 = t0, w = w,
              sv = sv, err_tol = eps, log = FALSE, switch_mech = "small",
              n_terms_small = "SWSE", summation_small = "2014"),
  dif = numeric(N),
  log_res = dfddm(rt = rt, response = "lower", a = a, v = v, t0 = t0, w = w,
                  sv = sv, err_tol = eps, log = TRUE, switch_mech = "small",
                  n_terms_small = "SWSE", summation_small = "2014")
)
SWSE_t_17 <- data.frame(
  res = dfddm(rt = rt, response = "lower", a = a, v = v, t0 = t0, w = w,
              sv = sv, err_tol = eps, log = FALSE, switch_mech = "eff_rt",
              n_terms_small = "SWSE", summation_small = "2017"),
  dif = numeric(N),
  log_res = dfddm(rt = rt, response = "lower", a = a, v = v, t0 = t0, w = w,
                  sv = sv, err_tol = eps, log = TRUE, switch_mech = "eff_rt",
                  n_terms_small = "SWSE", summation_small = "2017")
)
SWSE_t_14 <- data.frame(
  res = dfddm(rt = rt, response = "lower", a = a, v = v, t0 = t0, w = w,
              sv = sv, err_tol = eps, log = FALSE, switch_mech = "eff_rt",
              n_terms_small = "SWSE", summation_small = "2014"),
  dif = numeric(N),
  log_res = dfddm(rt = rt, response = "lower", a = a, v = v, t0 = t0, w = w,
                  sv = sv, err_tol = eps, log = TRUE, switch_mech = "eff_rt",
                  n_terms_small = "SWSE", summation_small = "2014")
)
SWSE_b_17 <- data.frame(
  res = dfddm(rt = rt, response = "lower", a = a, v = v, t0 = t0, w = w,
              sv = sv, err_tol = eps, log = FALSE, switch_mech = "terms_large",
              n_terms_small = "SWSE", summation_small = "2017"),
  dif = numeric(N),
  log_res = dfddm(rt = rt, response = "lower", a = a, v = v, t0 = t0, w = w,
                  sv = sv, err_tol = eps, log = TRUE,
                  switch_mech = "terms_large", n_terms_small = "SWSE",
                  summation_small = "2017")
)
SWSE_b_14 <- data.frame(
  res = dfddm(rt = rt, response = "lower", a = a, v = v, t0 = t0, w = w,
              sv = sv, err_tol = eps, log = FALSE, switch_mech = "terms_large",
              n_terms_small = "SWSE", summation_small = "2014"),
  dif = numeric(N),
  log_res = dfddm(rt = rt, response = "lower", a = a, v = v, t0 = t0, w = w,
                  sv = sv, err_tol = eps, log = TRUE,
                  switch_mech = "terms_large", n_terms_small = "SWSE",
                  summation_small = "2014")
)
Gondan_s_17 <- data.frame(
  res = dfddm(rt = rt, response = "lower", a = a, v = v, t0 = t0, w = w,
              sv = sv, err_tol = eps, log = FALSE, switch_mech = "small",
              n_terms_small = "Gondan", summation_small = "2017"),
  dif = numeric(N),
  log_res = dfddm(rt = rt, response = "lower", a = a, v = v, t0 = t0, w = w,
                  sv = sv, err_tol = eps, log = TRUE, switch_mech = "small",
                  n_terms_small = "Gondan", summation_small = "2017")
)
Gondan_s_14 <- data.frame(
  res = dfddm(rt = rt, response = "lower", a = a, v = v, t0 = t0, w = w,
              sv = sv, err_tol = eps, log = FALSE, switch_mech = "small",
              n_terms_small = "Gondan", summation_small = "2014"),
  dif = numeric(N),
  log_res = dfddm(rt = rt, response = "lower", a = a, v = v, t0 = t0, w = w,
                  sv = sv, err_tol = eps, log = TRUE, switch_mech = "small",
                  n_terms_small = "Gondan", summation_small = "2014")
)
Gondan_b_17 <- data.frame(
  res = dfddm(rt = rt, response = "lower", a = a, v = v, t0 = t0, w = w,
              sv = sv, err_tol = eps, log = FALSE, switch_mech = "terms",
              n_terms_small = "Gondan", summation_small = "2017"),
  dif = numeric(N),
  log_res = dfddm(rt = rt, response = "lower", a = a, v = v, t0 = t0, w = w,
                  sv = sv, err_tol = eps, log = TRUE, switch_mech = "terms",
                  n_terms_small = "Gondan", summation_small = "2017")
)
Gondan_b_14 <- data.frame(
  res = dfddm(rt = rt, response = "lower", a = a, v = v, t0 = t0, w = w,
              sv = sv, err_tol = eps, log = FALSE, switch_mech = "terms",
              n_terms_small = "Gondan", summation_small = "2014"),
  dif = numeric(N),
  log_res = dfddm(rt = rt, response = "lower", a = a, v = v, t0 = t0, w = w,
                  sv = sv, err_tol = eps, log = TRUE, switch_mech = "terms",
                  n_terms_small = "Gondan", summation_small = "2014")
)
Navarro_s_17 <- data.frame(
  res = dfddm(rt = rt, response = "lower", a = a, v = v, t0 = t0, w = w,
              sv = sv, err_tol = eps, log = FALSE, switch_mech = "small",
              n_terms_small = "Navarro", summation_small = "2017"),
  dif = numeric(N),
  log_res = dfddm(rt = rt, response = "lower", a = a, v = v, t0 = t0, w = w,
                  sv = sv, err_tol = eps, log = TRUE, switch_mech = "small",
                  n_terms_small = "Navarro", summation_small = "2017")
)
Navarro_s_14 <- data.frame(
  res = dfddm(rt = rt, response = "lower", a = a, v = v, t0 = t0, w = w,
              sv = sv, err_tol = eps, log = FALSE, switch_mech = "small",
              n_terms_small = "Navarro", summation_small = "2014"),
  dif = numeric(N),
  log_res = dfddm(rt = rt, response = "lower", a = a, v = v, t0 = t0, w = w,
                  sv = sv, err_tol = eps, log = TRUE, switch_mech = "small",
                  n_terms_small = "Navarro", summation_small = "2014")
)
Navarro_b_17 <- data.frame(
  res = dfddm(rt = rt, response = "lower", a = a, v = v, t0 = t0, w = w,
              sv = sv, err_tol = eps, log = FALSE, switch_mech = "terms",
              n_terms_small = "Navarro", summation_small = "2017"),
  dif = numeric(N),
  log_res = dfddm(rt = rt, response = "lower", a = a, v = v, t0 = t0, w = w,
                  sv = sv, err_tol = eps, log = TRUE, switch_mech = "terms",
                  n_terms_small = "Navarro", summation_small = "2017")
)
Navarro_b_14 <- data.frame(
  res = dfddm(rt = rt, response = "lower", a = a, v = v, t0 = t0, w = w,
              sv = sv, err_tol = eps, log = FALSE, switch_mech = "terms",
              n_terms_small = "Navarro", summation_small = "2014"),
  dif = numeric(N),
  log_res = dfddm(rt = rt, response = "lower", a = a, v = v, t0 = t0, w = w,
                  sv = sv, err_tol = eps, log = TRUE, switch_mech = "terms",
                  n_terms_small = "Navarro", summation_small = "2014")
)
Navarro_l <- data.frame(
  res = dfddm(rt = rt, response = "lower", a = a, v = v, t0 = t0, w = w,
              sv = sv, err_tol = eps, log = FALSE, switch_mech = "large",
              n_terms_small = "Navarro"),
  dif = numeric(N),
  log_res = dfddm(rt = rt, response = "lower", a = a, v = v, t0 = t0, w = w,
                  sv = sv, err_tol = eps, log = TRUE, switch_mech = "large",
                  n_terms_small = "Navarro")
)

# non-fddm methods
t <- rt - t0
M <- exp(v * a * w + v*v * t / 2 +
         (sv*sv * a*a * w*w - 2 * v * a * w - v*v * t) / (2 + 2 * sv*sv * t)
        ) / sqrt(1 + sv*sv * t)
if (require("RWiener")) {
  RWiener  <- data.frame(
    res = numeric(N),
    dif = numeric(N)
  )
  for (i in 1:N) { # RWiener can't handle model parameters as vectors
    RWiener[i, "res"] <- dwiener(rt[i], resp = "lower", alpha = a[i],
                                 delta = v[i], tau = t0, beta = w[i],
                                 give_log = FALSE)
  }
  RWiener[["res"]] <- M * RWiener[["res"]]
}
Gondan_R <- data.frame(
  res = M * fs(t = t, a = a, v = v, w = w, eps = eps),
  dif = numeric(N)
)
if (require("rtdists")) {
  rtdists <- data.frame(
    res = ddiffusion(rt, "lower", a = a, v = v, t0 = t0, z = w*a, sv = sv),
    dif = numeric(N)
  )
}



# Calculate differences (use fddm's SWSE_t_17 method as truth)
ans <- SWSE_t_17[["res"]]
SWSE_s_17[["dif"]] <- abs(SWSE_s_17[["res"]] - ans)
SWSE_s_14[["dif"]] <- abs(SWSE_s_14[["res"]] - ans)
SWSE_t_17[["dif"]] <- abs(SWSE_t_17[["res"]] - ans)
SWSE_t_14[["dif"]] <- abs(SWSE_t_14[["res"]] - ans)
SWSE_b_17[["dif"]] <- abs(SWSE_b_17[["res"]] - ans)
SWSE_b_14[["dif"]] <- abs(SWSE_b_14[["res"]] - ans)
Gondan_s_17[["dif"]] <- abs(Gondan_s_17[["res"]] - ans)
Gondan_s_14[["dif"]] <- abs(Gondan_s_14[["res"]] - ans)
Gondan_b_17[["dif"]] <- abs(Gondan_b_17[["res"]] - ans)
Gondan_b_14[["dif"]] <- abs(Gondan_b_14[["res"]] - ans)
Navarro_s_17[["dif"]] <- abs(Navarro_s_17[["res"]] - ans)
Navarro_s_14[["dif"]] <- abs(Navarro_s_14[["res"]] - ans)
Navarro_b_17[["dif"]] <- abs(Navarro_b_17[["res"]] - ans)
Navarro_b_14[["dif"]] <- abs(Navarro_b_14[["res"]] - ans)
Navarro_l[["dif"]] <- abs(Navarro_l[["res"]] - ans)
if (require("RWiener")) {
  RWiener[["dif"]] <- abs(RWiener[["res"]] - ans)
}
Gondan_R[["dif"]] <- abs(Gondan_R[["res"]] - ans)
if (require("rtdists")) {
  rtdists[["dif"]] <- abs(rtdists[["res"]] - ans)
}





#---------------------- Testing -----------------------------------------------#
# Ensure all densities are non-negative
test_that("Non-negativity of densities", {
  expect_true(all(SWSE_s_17[["res"]] >= 0))
  expect_true(all(SWSE_s_14[["res"]] >= 0))
  expect_true(all(SWSE_t_17[["res"]] >= 0))
  expect_true(all(SWSE_t_14[["res"]] >= 0))
  expect_true(all(SWSE_b_17[["res"]] >= 0))
  expect_true(all(SWSE_b_14[["res"]] >= 0))
  expect_true(all(Gondan_s_17[["res"]] >= 0))
  expect_true(all(Gondan_s_14[["res"]] >= 0))
  expect_true(all(Gondan_b_17[["res"]] >= 0))
  expect_true(all(Gondan_b_14[["res"]] >= 0))
  expect_true(all(Navarro_s_17[["res"]] >= 0))
  expect_true(all(Navarro_s_14[["res"]] >= 0))
  expect_true(all(Navarro_b_17[["res"]] >= 0))
  expect_true(all(Navarro_b_14[["res"]] >= 0))
  expect_true(all(Navarro_l[["res"]] >= 0))
  if (require("RWiener")) {
    expect_true(all(RWiener[["res"]] >= 0))
  }
  expect_true(all(Gondan_R[["res"]] >= -eps)) # density between 0 and -eps := 0
  if (require("rtdists")) {
    expect_true(all(rtdists[["res"]] >= 0))
  }
})



# Test accuracy within 2*eps (allows for convergence from above and below)
test_that("Consistency among internal methods", {
  expect_true(all(SWSE_s_17[["dif"]] <= 2 * eps))
  expect_true(all(SWSE_s_14[["dif"]] <= 2 * eps))
  expect_true(all(SWSE_t_17[["dif"]] <= 2 * eps))
  expect_true(all(SWSE_t_14[["dif"]] <= 2 * eps))
  expect_true(all(SWSE_b_17[["dif"]] <= 2 * eps))
  expect_true(all(SWSE_b_14[["dif"]] <= 2 * eps))
  expect_true(all(Gondan_s_17[["dif"]] <= 2 * eps))
  expect_true(all(Gondan_s_14[["dif"]] <= 2 * eps))
  expect_true(all(Gondan_b_17[["dif"]] <= 2 * eps))
  expect_true(all(Gondan_b_14[["dif"]] <= 2 * eps))
  expect_true(all(Navarro_s_17[["dif"]] <= 2 * eps))
  expect_true(all(Navarro_s_14[["dif"]] <= 2 * eps))
  expect_true(all(Navarro_b_17[["dif"]] <= 2 * eps))
  expect_true(all(Navarro_b_14[["dif"]] <= 2 * eps))
  testthat::skip_on_os("solaris")
  testthat::skip_if(dfddm(rt = 0.001, response = "lower",
                          a = 5, v = -5, t0 = 1e-4, w = 0.8, sv = 1.5,
                          err_tol = 1e-6, log = FALSE, switch_mech = "large") >
                    1e-6)
  expect_true(all(Navarro_l[rt/a/a >= 0.009, "dif"] < 2 * eps)) # see KE 1
})

test_that("Accuracy relative to established packages", {
  if (require("RWiener")) {
    expect_true(all(RWiener[sv < SV_THRESH, "dif"] <= 2 * eps)) # see KE 2
  }
  if (require("rtdists")) {
    expect_true(all(rtdists[["dif"]] <= 2 * eps))
  }
  testthat::skip_on_os("solaris")
  testthat::skip_if(dfddm(rt = 0.001, response = "lower",
                          a = 5, v = -5, t0 = 1e-4, w = 0.8, sv = 1.5,
                          err_tol = 1e-6, log = FALSE, switch_mech = "large") >
                    1e-6)
  expect_true(all(Gondan_R[sv < SV_THRESH, "dif"] <= 2 * eps)) # see KE 2
})



# Test consistency in fddm log vs non-log (see KE 3)
test_that("Log-Consistency among internal methods", {
  expect_equal(SWSE_s_17[SWSE_s_17[["res"]] > eps*eps, "log_res"],
               log(SWSE_s_17[SWSE_s_17[["res"]] > eps*eps, "res"]))
  expect_equal(SWSE_s_14[SWSE_s_14[["res"]] > eps*eps, "log_res"],
               log(SWSE_s_14[SWSE_s_14[["res"]] > eps*eps, "res"]))
  expect_equal(SWSE_t_17[SWSE_t_17[["res"]] > eps*eps, "log_res"],
               log(SWSE_t_17[SWSE_t_17[["res"]] > eps*eps, "res"]))
  expect_equal(SWSE_t_14[SWSE_t_14[["res"]] > eps*eps, "log_res"],
               log(SWSE_t_14[SWSE_t_14[["res"]] > eps*eps, "res"]))
  expect_equal(SWSE_b_17[SWSE_b_17[["res"]] > eps*eps, "log_res"],
               log(SWSE_b_17[SWSE_b_17[["res"]] > eps*eps, "res"]))
  expect_equal(SWSE_b_14[SWSE_b_14[["res"]] > eps*eps, "log_res"],
               log(SWSE_b_14[SWSE_b_14[["res"]] > eps*eps, "res"]))
  expect_equal(Gondan_s_17[Gondan_s_17[["res"]] > eps*eps, "log_res"],
               log(Gondan_s_17[Gondan_s_17[["res"]] > eps*eps, "res"]))
  expect_equal(Gondan_s_14[Gondan_s_14[["res"]] > eps*eps, "log_res"],
               log(Gondan_s_14[Gondan_s_14[["res"]] > eps*eps, "res"]))
  expect_equal(Gondan_b_17[Gondan_b_17[["res"]] > eps*eps, "log_res"],
               log(Gondan_b_17[Gondan_b_17[["res"]] > eps*eps, "res"]))
  expect_equal(Gondan_b_14[Gondan_b_14[["res"]] > eps*eps, "log_res"],
               log(Gondan_b_14[Gondan_b_14[["res"]] > eps*eps, "res"]))
  expect_equal(Navarro_s_17[Navarro_s_17[["res"]] > eps*eps, "log_res"],
               log(Navarro_s_17[Navarro_s_17[["res"]] > eps*eps, "res"]))
  expect_equal(Navarro_s_14[Navarro_s_14[["res"]] > eps*eps, "log_res"],
               log(Navarro_s_14[Navarro_s_14[["res"]] > eps*eps, "res"]))
  expect_equal(Navarro_b_17[Navarro_b_17[["res"]] > eps*eps, "log_res"],
               log(Navarro_b_17[Navarro_b_17[["res"]] > eps*eps, "res"]))
  expect_equal(Navarro_b_14[Navarro_b_14[["res"]] > eps*eps, "log_res"],
               log(Navarro_b_14[Navarro_b_14[["res"]] > eps*eps, "res"]))
  expect_equal(Navarro_l[Navarro_l[["res"]] > eps*eps, "log_res"],
               log(Navarro_l[Navarro_l[["res"]] > eps*eps, "res"]))
})





### Known Errors (KE) ###
#
# 1) The "large-time" variant is unstable for small effective response times
#    ( (rt - t0) / (a*a) < 0.009 ) and produces inaccurate densities.
#
# 2) Both RWiener and Gondan_R divide the error tolerance by the multiplicative
#    term outside of the summation. Since the outside term is different when
#    $sv > 0$, the approximations use the incorrect error tolerance for
#    $sv > 0$. This affects the number of terms required in the summation to
#    achieve the desired precision, thus not actually achieving that desired
#    precision. This issue is fixed in our implementation of the Gondan method,
#    `switch_mech = "small"`, `n_terms_small = "Gondan"`. For an example of this
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
# 3) When calculating the log of the density, it is better to use the built-in
#    log option. For very small densities, simply calculating the density can
#    cause rounding issues that result in a density of zero (thus the log of the
#    density becomes -Inf). Using the built-in log option avoids some of these
#    rounding issues by exploiting the algebraic properties of the logarithm.
#    Also note that sometimes the densities are just too small (i.e. extremely
#    negative) and the logarithm function returns a value of -Inf, so we discard
#    the samples whose density is very small (less than eps*eps = 1e-12).
