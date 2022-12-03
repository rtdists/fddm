context("Testing the partial derivatives of the density function for accuracy")
# note: allowance is 4*err because WienR doesn't account for 2 sums

test_that("Accuracy relative to WienR", {
  testthat::skip_if_not_installed("WienR")
  library("WienR")

  ### Define Parameter Space ###
  RT <- c(0.001, 0.1, 1, 2, 3, 4, 5, 10, 30)
  A <- c(0.25, 0.5, 1, 2.5, 5)
  V <- c(-5, -2, 0, 2, 5)
  W <- c(0.2, 0.5, 0.8)
  SV <- c(0, 0.5, 1, 1.5)
  t0 <- 0
  err <- 1e-12 # default setting from WienR

  nRT <- length(RT)
  nA <- length(A)
  nV <- length(V)
  nW <- length(W)
  nSV <- length(SV)
  N <- nRT * nA * nV * nW * nSV

  resp <- "lower"
  respn <- rep(resp, N)
  rt <- rep(RT, each = nSV * nW * nV * nA, times = 1)
  a <-  rep(A,  each = nSV * nW * nV, times = nRT)
  v <-  rep(V,  each = nSV * nW, times = nRT * nA)
  w <-  rep(W,  each = nSV, times = nRT * nA * nV)
  sv <- rep(SV, each = 1, times = nRT * nA * nV * nW)



  ###
  ### dv ###
  df <- data.frame(
    rt = rt - t0,
    a = a,
    v = v,
    w = w,
    sv = sv,
    res_WienR = dvWienerPDF(t = rt, response = respn, a = a, v = v, w = w,
                            t0 = t0, sv = sv)[["deriv"]],
    res_fddm = dv_dfddm(rt = rt, response = resp, a = a, v = v, t0 = t0,
                        w = w, sv = sv, err_tol = err),
    diff = numeric(N)
  )
  df[["diff"]] <- df[["res_WienR"]] - df[["res_fddm"]]

  # agreement with WienR
  expect_true(all(abs(df[["diff"]]) <= 4*err))



  ###
  ### da ###
  df[["res_WienR"]] <- daWienerPDF(t = rt, response = respn, a = a, v = v,
                                   w = w, t0 = t0, sv = sv)[["deriv"]]
  df[["res_fddm"]] <- da_dfddm(rt = rt, response = resp, a = a, v = v, t0 = t0,
                               w = w, sv = sv, err_tol = err)
  df[["diff"]] <- df[["res_WienR"]] - df[["res_fddm"]]

  # agreement with WienR
  expect_true(sum(abs(df[["diff"]]) <= 4*err) / nrow(df) > 0.98)

  # disagreements with WienR
  df_bad <- df[abs(df[["diff"]]) > 4*err, ]

  # small-time and large-time in fddm are equal
  fddm_sm <- da_dfddm(rt = df_bad[["rt"]], response = resp,
                      a = df_bad[["a"]], v = df_bad[["v"]], w = df_bad[["w"]],
                      sv = df_bad[["sv"]], t0 = 0, sl_thresh = Inf,
                      err_tol = err)
  fddm_lg <- da_dfddm(rt = df_bad[["rt"]], response = resp,
                      a = df_bad[["a"]], v = df_bad[["v"]], w = df_bad[["w"]],
                      sv = df_bad[["sv"]], t0 = 0, sl_thresh = -1,
                      err_tol = err)
  expect_true(sum(abs(fddm_sm - fddm_lg) <= 2*err) / nrow(df_bad) == 1)

  # disagreements are small relative to the calculation
  # WienR precision not guaranteed for sv > 0
  df_bad <- df_bad[df_bad[["sv"]] == 0, ]
  expect_true(all(abs(df_bad[["diff"]]) / abs(df_bad[["res_WienR"]]) <=
                  err*4e-2))



  ###
  ### dt ###
  df[["res_WienR"]] <- dtWienerPDF(t = rt, response = respn, a = a, v = v,
                                   w = w, t0 = t0, sv = sv)[["deriv"]]
  df[["res_fddm"]] <- dt_dfddm(rt = rt, response = resp, a = a, v = v, t0 = t0,
                               w = w, sv = sv, err_tol = err)
  df[["diff"]] <- df[["res_WienR"]] - df[["res_fddm"]]

  # agreement with WienR
  expect_true(sum(abs(df[["diff"]]) <= 4*err) / nrow(df) > 0.97)

  # disagreements with WienR (WienR precision not guaranteed for sv > 0)
  df_bad <- df[abs(df[["diff"]]) > 4*err & df[["sv"]] == 0, ]

  # disagreements are small relative to the calculation
  expect_true(all(abs(df_bad[["diff"]]) / abs(df_bad[["res_WienR"]]) <=
                  err*4e-2))



  ###
  ### dt0 ###
  df[["res_WienR"]] <- dt0WienerPDF(t = rt, response = respn, a = a, v = v,
                                    w = w, t0 = t0, sv = sv)[["deriv"]]
  df[["res_fddm"]] <- dt0_dfddm(rt = rt, response = resp, a = a, v = v, t0 = t0,
                                w = w, sv = sv, err_tol = err)
  df[["diff"]] <- df[["res_WienR"]] - df[["res_fddm"]]

  # agreement with WienR
  expect_true(sum(abs(df[["diff"]]) <= 4*err) / nrow(df) > 0.97)

  # disagreements with WienR (WienR precision not guaranteed for sv > 0)
  df_bad <- df[abs(df[["diff"]]) > 4*err & df[["sv"]] == 0, ]

  # disagreements are small relative to the calculation
  expect_true(all(abs(df_bad[["diff"]]) / abs(df_bad[["res_WienR"]]) <=
                  err*4e-2))



  ###
  ### dw ###
  df[["res_WienR"]] <- dwWienerPDF(t = rt, response = respn, a = a, v = v,
                                   w = w, t0 = t0, sv = sv)[["deriv"]]
  df[["res_fddm"]] <- dw_dfddm(rt = rt, response = resp, a = a, v = v, t0 = t0,
                               w = w, sv = sv, err_tol = err)
  df[["diff"]] <- df[["res_WienR"]] - df[["res_fddm"]]

  # agreement with WienR
  expect_true(sum(abs(df[["diff"]]) <= 4*err) / nrow(df) > 0.98)

  # disagreements with WienR
  df_bad <- df[abs(df[["diff"]]) > 4*err, ]

  # small-time and large-time in fddm are equal
  fddm_sm <- dw_dfddm(rt = df_bad[["rt"]], response = resp,
                      a = df_bad[["a"]], v = df_bad[["v"]], w = df_bad[["w"]],
                      sv = df_bad[["sv"]], t0 = 0, sl_thresh = -1,
                      err_tol = err)
  fddm_lg <- dw_dfddm(rt = df_bad[["rt"]], response = resp,
                      a = df_bad[["a"]], v = df_bad[["v"]], w = df_bad[["w"]],
                      sv = df_bad[["sv"]], t0 = 0, sl_thresh = Inf,
                      err_tol = err)
  expect_true(sum(abs(fddm_sm - fddm_lg) <= 2*err) / nrow(df_bad) == 1)



  ###
  ### dsv ###
  SV <- c(0.5, 1, 1.5) # remove sv = 0
  nSV <- length(SV)
  N <- nRT * nA * nV * nW * nSV
  rt <- rep(RT, each = nSV * nW * nV * nA, times = 1) - t0
  a <-  rep(A,  each = nSV * nW * nV, times = nRT)
  v <-  rep(V,  each = nSV * nW, times = nRT * nA)
  w <-  rep(W,  each = nSV, times = nRT * nA * nV)
  sv <- rep(SV, each = 1, times = nRT * nA * nV * nW)
  N <- nRT * nA * nV * nW * nSV
  respn <- rep(resp, N)

  df <- data.frame(
    rt = rt - t0,
    a = a,
    v = v,
    w = w,
    sv = sv,
    res_WienR = dsvWienerPDF(t = rt, response = respn, a = a, v = v,
                             w = w, t0 = t0, sv = sv)[["deriv"]],
    res_fddm = dsv_dfddm(rt = rt, response = resp, a = a, v = v, t0 = t0,
                         w = w, sv = sv, err_tol = err),
    diff = numeric(N)
  )
  df[["diff"]] <- df[["res_WienR"]] - df[["res_fddm"]]

  # agreement with WienR
  expect_true(all(abs(df[["diff"]]) <= 4*err))
})





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
  err <- 1e-12
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
  ### dv ###
  df <- data.frame(
    rt = rt,
    a = a,
    rtaa = rt/(a*a),
    v = v,
    w = w,
    sv = sv,
    res_small = dv_dfddm(rt = rt, response = resp, v = v, a = a, t0 = t0,
                         w = w, sv = sv, sl_thresh = 1000, err_tol = err),
    res_large = dv_dfddm(rt = rt, response = resp, v = v, a = a, t0 = t0,
                         w = w, sv = sv, sl_thresh = 0, err_tol = err),
    res_appx = numeric(N),
    diff_small = numeric(N),
    diff_large = numeric(N)
  )
  dv_wrap <- function(v, p) {
    return(dfddm(rt = p[1], response = resp, v = v, a = p[2], t0 = 0,
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
  ### da ###
  df[["res_small"]] <- da_dfddm(rt = rt, response = resp, v = v, a = a, t0 = t0,
                                w = w, sv = sv, sl_thresh = 1000, err_tol = err)
  df[["res_large"]] <- da_dfddm(rt = rt, response = resp, v = v, a = a, t0 = t0,
                                w = w, sv = sv, sl_thresh = 0, err_tol = err)
  da_wrap <- function(a, p) {
    return(dfddm(rt = p[1], response = resp, v = p[2], a = a, t0 = 0,
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
  expect_true(all(abs(df[df[["rtaa"]] > 4e-03, "diff_large"]) <= eps))



  ###
  ### dt ###
  df[["res_small"]] <- dt_dfddm(rt = rt, response = resp, v = v, a = a, t0 = t0,
                                w = w, sv = sv, sl_thresh = 1000, err_tol = err)
  df[["res_large"]] <- dt_dfddm(rt = rt, response = resp, v = v, a = a, t0 = t0,
                                w = w, sv = sv, sl_thresh = 0, err_tol = err)
  dt_wrap <- function(rt, p) {
    return(dfddm(rt = rt, response = resp, v = p[1], a = p[2], t0 = 0,
                 w = p[3], sv = p[4], err_tol = err))
  }
  for (i in seq_len(N)) {
    df[i, "res_appx"] <-
      numDeriv::grad(func = dt_wrap, x = c("rt" = rt[i]), method = "Richardson",
                     p = c(v = v[i], a = a[i], w = w[i], sv = sv[i]))
  }
  df[["diff_small"]] <- df[["res_appx"]] - df[["res_small"]]
  df[["diff_large"]] <- df[["res_appx"]] - df[["res_large"]]

  expect_true(sum(abs(df[["diff_small"]]) <= eps) / N > 0.99)
  expect_true(sum(abs(df[["diff_large"]]) <= eps) / N > 0.98)

  # disagreements are small relative to the calculation
  expect_true(all(abs(df[abs(df[["diff_small"]]) > eps, "diff_small"]) /
                  abs(df[abs(df[["diff_small"]]) > eps, "res_appx"]) <= eps))
  expect_true(all(abs(df[abs(df[["diff_large"]]) > eps & df[["rtaa"]] > 0.01,
                         "diff_small"]) /
                  abs(df[abs(df[["diff_large"]]) > eps & df[["rtaa"]] > 0.01,
                         "res_appx"]) <= eps))



  ###
  ### dt0 ###
  df[["res_small"]] <- dt0_dfddm(rt = rt, response = resp, v = v, a = a,
                                 t0 = t0, w = w, sv = sv, sl_thresh = 1000,
                                 err_tol = err)
  df[["res_large"]] <- dt0_dfddm(rt = rt, response = resp, v = v, a = a,
                                 t0 = t0, w = w, sv = sv, sl_thresh = 0,
                                 err_tol = err)
  dt0_wrap <- function(rt, p) {
    return(dfddm(rt = rt, response = resp, v = p[1], a = p[2], t0 = 0,
                 w = p[3], sv = p[4], err_tol = err))
  }
  for (i in seq_len(N)) {
    df[i, "res_appx"] <- -1 *
      numDeriv::grad(func = dt0_wrap, x = c("rt" = rt[i]),
                     method = "Richardson",
                     p = c(v = v[i], a = a[i], w = w[i], sv = sv[i]))
  }
  df[["diff_small"]] <- df[["res_appx"]] - df[["res_small"]]
  df[["diff_large"]] <- df[["res_appx"]] - df[["res_large"]]

  expect_true(sum(abs(df[["diff_small"]]) <= eps) / N > 0.99)
  expect_true(sum(abs(df[["diff_large"]]) <= eps) / N > 0.98)

  # disagreements are small relative to the calculation
  expect_true(all(abs(df[abs(df[["diff_small"]]) > eps, "diff_small"]) /
                  abs(df[abs(df[["diff_small"]]) > eps, "res_appx"]) <= eps))
  expect_true(all(abs(df[abs(df[["diff_large"]]) > eps & df[["rtaa"]] > 0.01,
                         "diff_small"]) /
                  abs(df[abs(df[["diff_large"]]) > eps & df[["rtaa"]] > 0.01,
                         "res_appx"]) <= eps))



  ###
  ### dw ###
  df[["res_small"]] <- dw_dfddm(rt = rt, response = resp, v = v, a = a, t0 = t0,
                                w = w, sv = sv, sl_thresh = 1000, err_tol = err)
  df[["res_large"]] <- dw_dfddm(rt = rt, response = resp, v = v, a = a, t0 = t0,
                                w = w, sv = sv, sl_thresh = 0, err_tol = err)
  dw_wrap <- function(w, p) {
    return(dfddm(rt = p[1], response = resp, v = p[2], a = p[3], t0 = 0,
                 w = w, sv = p[4], err_tol = err))
  }
  for (i in seq_len(N)) {
    df[i, "res_appx"] <-
      numDeriv::grad(func = dw_wrap, x = c("w" = w[i]), method = "Richardson",
                    p = c(rt = rt[i], v = v[i], a = a[i], sv = sv[i]))
  }
  df[["diff_small"]] <- df[["res_appx"]] - df[["res_small"]]
  df[["diff_large"]] <- df[["res_appx"]] - df[["res_large"]]

  expect_true(all(abs(df[["diff_small"]]) <= eps))
  expect_true(all(abs(df[["diff_large"]]) <= eps))



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
    return(dfddm(rt = p[1], response = resp, v = p[2], a = p[3], t0 = 0,
                 w = p[4], sv = sv, err_tol = err))
  }
  df <- data.frame(
    rt = rt,
    a = a,
    v = v,
    w = w,
    sv = sv,
    res_small = dsv_dfddm(rt = rt, response = resp, v = v, a = a, t0 = t0,
                          w = w, sv = sv, sl_thresh = 1000, err_tol = err),
    res_large = dsv_dfddm(rt = rt, response = resp, v = v, a = a, t0 = t0,
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
