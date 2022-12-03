context("Testing the second order partial derivatives of the density function
for accuracy")
# note: eps = 2*sqrt(err) because we're testing against numerical approximations



test_that("Accuracy relative to numerical approximations", {
  testthat::skip_if_not_installed("numDeriv")
  library("numDeriv")

  ### Define Parameter Space ###
  RT <- c(0.001, 0.1, 1, 2, 3, 4, 5, 10, 30)
  A <- c(0.25, 0.5, 1, 2.5, 5)
  V <- c(-5, -2, 0, 2, 5)
  W <- c(0.2, 0.5, 0.8)
  SV <- c(0, 0.5, 1, 1.5)
  t0 <- 0
  err <- 1e-6
  eps <- 2 * sqrt(err)

  nRT <- length(RT)
  nA <- length(A)
  nV <- length(V)
  nW <- length(W)
  nSV <- length(SV)
  N <- nRT * nA * nV * nW * nSV

  resp <- "lower"
  rt <- rep(RT, each = nSV * nW * nV * nA, times = 1) - t0
  a <-  rep(A,  each = nSV * nW * nV, times = nRT)
  v <-  rep(V,  each = nSV * nW, times = nRT * nA)
  w <-  rep(W,  each = nSV, times = nRT * nA * nV)
  sv <- rep(SV, each = 1, times = nRT * nA * nV * nW)



  ###
  ### dv2 ###
  df <- data.frame(
    rt = rt,
    a = a,
    rtaa = rt/(a*a),
    v = v,
    w = w,
    sv = sv,
    res_small = dv2_dfddm(rt = rt, response = resp, v = v, a = a, t0 = t0,
                          w = w, sv = sv, sl_thresh = 1000, err_tol = err),
    res_large = dv2_dfddm(rt = rt, response = resp, v = v, a = a, t0 = t0,
                          w = w, sv = sv, sl_thresh = 0, err_tol = err),
    res_appx = numeric(N),
    diff_small = numeric(N),
    diff_large = numeric(N)
  )
  dv_wrap <- function(v, p) {
    return(dv_dfddm(rt = p[1], response = resp, v = v, a = p[2], t0 = 0,
                    w = p[3], sv = p[4], err_tol = err))
  }
  for (i in seq_len(N)) {
    df[i, "res_appx"] <-
      numDeriv::grad(func = dv_wrap, x = c("v" = v[i]), method = "Richardson",
                    p = c(rt = rt[i], a = a[i], w = w[i], sv = sv[i]))
  }
  df[["diff_small"]] <- df[["res_appx"]] - df[["res_small"]]
  df[["diff_large"]] <- df[["res_appx"]] - df[["res_large"]]

  expect_true(all(abs(df[["diff_small"]]) <= eps))
  expect_true(all(abs(df[["diff_large"]]) <= eps))



  ###
  ### da2 ###
  df[["res_small"]] <- da2_dfddm(rt = rt, response = resp, v = v, a = a,
                                 t0 = t0, w = w, sv = sv, sl_thresh = 1000,
                                 err_tol = err)
  df[["res_large"]] <- da2_dfddm(rt = rt, response = resp, v = v, a = a,
                                 t0 = t0, w = w, sv = sv, sl_thresh = 0,
                                 err_tol = err)
  da_wrap <- function(a, p) {
    return(da_dfddm(rt = p[1], response = resp, v = p[2], a = a, t0 = 0,
                    w = p[3], sv = p[4], err_tol = err))
  }
  for (i in seq_len(N)) {
    df[i, "res_appx"] <-
      numDeriv::grad(func = da_wrap, x = c("a" = a[i]), method = "Richardson",
                    p = c(rt = rt[i], v = v[i], w = w[i], sv = sv[i]))
  }
  df[["diff_small"]] <- df[["res_appx"]] - df[["res_small"]]
  df[["diff_large"]] <- df[["res_appx"]] - df[["res_large"]]

  expect_true(all(abs(df[["diff_small"]]) <= eps))
  expect_true(all(abs(df[df[["rtaa"]] > 0.009, "diff_large"]) <= eps))



  ###
  ### dt02 (dt = -dt0; dt^2 = dt0^2) ###
  df[["res_small"]] <- dt02_dfddm(rt = rt, response = resp, v = v, a = a,
                                  t0 = t0, w = w, sv = sv, sl_thresh = 1000,
                                  err_tol = err)
  df[["res_large"]] <- dt02_dfddm(rt = rt, response = resp, v = v, a = a,
                                  t0 = t0, w = w, sv = sv, sl_thresh = 0,
                                  err_tol = err)
  dt_wrap <- function(rt, p) {
    return(dt_dfddm(rt = rt, response = resp, v = p[1], a = p[2], t0 = 0,
                    w = p[3], sv = p[4], err_tol = err))
  }
  for (i in seq_len(N)) {
    df[i, "res_appx"] <-
      numDeriv::grad(func = dt_wrap, x = c("rt" = rt[i]), method = "Richardson",
                    p = c(v = v[i], a = a[i], w = w[i], sv = sv[i]))
  }
  df[["diff_small"]] <- df[["res_appx"]] - df[["res_small"]]
  df[["diff_large"]] <- df[["res_appx"]] - df[["res_large"]]

  expect_true(sum(abs(df[["diff_small"]]) <= eps) / nrow(df) > 0.99)
  expect_true(sum(abs(df[["diff_large"]]) <= eps) / nrow(df) > 0.8)

  # the relative error is ok (discrepancies only occur when result gets large)
  expect_true(all(abs(df[abs(df[["diff_small"]]) > eps, "diff_small"]) /
                  abs(df[abs(df[["diff_small"]]) > eps, "res_appx"]) <= eps))
  expect_true(all(abs(df[abs(df[["diff_large"]]) > eps & df[["rtaa"]] > 0.1,
                         "diff_small"]) /
                  abs(df[abs(df[["diff_large"]]) > eps & df[["rtaa"]] > 0.1,
                         "res_appx"]) <= eps))



  ###
  ### dw2 ###
  df[["res_small"]] <- dw2_dfddm(rt = rt, response = resp, v = v, a = a,
                                 t0 = t0, w = w, sv = sv, sl_thresh = 1000,
                                 err_tol = err)
  df[["res_large"]] <- dw2_dfddm(rt = rt, response = resp, v = v, a = a,
                                 t0 = t0, w = w, sv = sv, sl_thresh = 0,
                                 err_tol = err)
  dw_wrap <- function(w, p) {
    return(dw_dfddm(rt = p[1], response = resp, v = p[2], a = p[3], t0 = 0,
                    w = w, sv = p[4], err_tol = err))
  }
  for (i in seq_len(N)) {
    df[i, "res_appx"] <-
      numDeriv::grad(func = dw_wrap, x = c("w" = w[i]), method = "Richardson",
                    p = c(rt = rt[i], v = v[i], a = a[i], sv = sv[i]))
  }
  df[["diff_small"]] <- df[["res_appx"]] - df[["res_small"]]
  df[["diff_large"]] <- df[["res_appx"]] - df[["res_large"]]

  # t=1, a=2.5, v=2, w=0.8, sv=1 yields a weird spike in the numerical approx
  expect_equal(sum(abs(df[["diff_small"]]) > eps), 1)
  expect_equal(sum(abs(df[df[["rtaa"]] > 0.009, "diff_large"]) > eps), 1)



  ###
  ### dsv2 ###
  SV <- c(0.5, 1, 1.5) # remove sv = 0
  nSV <- length(SV)
  N <- nRT * nA * nV * nW * nSV
  rt <- rep(RT, each = nSV * nW * nV * nA, times = 1) - t0
  a <-  rep(A,  each = nSV * nW * nV, times = nRT)
  v <-  rep(V,  each = nSV * nW, times = nRT * nA)
  w <-  rep(W,  each = nSV, times = nRT * nA * nV)
  sv <- rep(SV, each = 1, times = nRT * nA * nV * nW)

  dsv_wrap <- function(sv, p) {
    return(dsv_dfddm(rt = p[1], response = resp, v = p[2], a = p[3], t0 = 0,
                     w = p[4], sv = sv, err_tol = err))
  }
  df <- data.frame(
    rt = rt,
    a = a,
    v = v,
    w = w,
    sv = sv,
    res_small = dsv2_dfddm(rt = rt, response = resp, v = v, a = a, t0 = t0,
                           w = w, sv = sv, sl_thresh = 1000, err_tol = err),
    res_large = dsv2_dfddm(rt = rt, response = resp, v = v, a = a, t0 = t0,
                           w = w, sv = sv, sl_thresh = 0, err_tol = err),
    res_appx = numeric(N),
    diff_small = numeric(N),
    diff_large = numeric(N)
  )
  for (i in seq_len(N)) {
    df[i, "res_appx"] <-
      numDeriv::grad(func = dsv_wrap, x = c("sv" = sv[i]),
                     method = "Richardson",
                     p = c(rt = rt[i], v = v[i], a = a[i], w = w[i]))
  }
  df[["diff_small"]] <- df[["res_appx"]] - df[["res_small"]]
  df[["diff_large"]] <- df[["res_appx"]] - df[["res_large"]]

  expect_true(all(abs(df[["diff_small"]]) <= eps))
  expect_true(all(abs(df[["diff_large"]]) <= eps))
})
