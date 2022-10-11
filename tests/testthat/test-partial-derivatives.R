context("Testing the partial derivatives of the density function for accuracy")
# note: allowance is 4*err because WienR doesn't account for 2 sums

test_that("Accuracy relative to WienR", {
  testthat::skip_if_not_installed("WienR")
  library("WienR")

  ### Define Parameter Space ###
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
  err <- 1e-12 # default setting from WienR

  nRT <- length(RT)
  nA <- length(A)
  nV <- length(V)
  nW <- length(W)
  nSV <- length(SV)
  N <- nRT * nA * nV * nW * nSV

  resp <- rep("lower", N)
  rt <- rep(RT, each = nSV * nW * nV * nA, times = 1)
  a <-  rep(A,  each = nSV * nW * nV, times = nRT)
  v <-  rep(V,  each = nSV * nW, times = nRT * nA)
  w <-  rep(W,  each = nSV, times = nRT * nA * nV)
  sv <- rep(SV, each = 1, times = nRT * nA * nV * nW)



  ### dt ###
  df <- data.frame(
    rt = rt - t0,
    a = a,
    v = v,
    w = w,
    sv = sv,
    res_WienR = dtWienerPDF(t = rt, response = resp, a = a, v = v, w = w,
                            t0 = t0, sv = sv)[["deriv"]],
    res_fddm = dt_dfddm(rt = rt, response = resp, a = a, v = v, t0 = t0, w = w,
                        sv = sv, err_tol = err),
    diff = numeric(N)
  )
  df[["diff"]] <- df[["res_WienR"]] - df[["res_fddm"]]

  # agreement with WienR
  expect_true(sum(abs(df[["diff"]]) <= 4*err) / nrow(df) > 0.97)

  # disagreements with WienR (WienR precision not guaranteed for sv > 0)
  df_bad <- df[abs(df[["diff"]]) > 4*err & df[["sv"]] == 0, ]

  # disagreements are small relative to the calculation
  expect_true(all(abs(df_bad[["diff"]]) / abs(df_bad[["res_WienR"]]) <=
                  err*4e-2))



  ### dt0 ###
  df[["res_WienR"]] <- dt0WienerPDF(t = rt, response = resp, a = a, v = v,
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



  ### da ###
  df[["res_WienR"]] <- daWienerPDF(t = rt, response = resp, a = a, v = v, w = w,
                                  t0 = t0, sv = sv)[["deriv"]]
  df[["res_fddm"]] <- da_dfddm(rt = rt, response = resp, a = a, v = v, t0 = t0,
                              w = w, sv = sv, err_tol = err)
  df[["diff"]] <- df[["res_WienR"]] - df[["res_fddm"]]

  # agreement with WienR
  expect_true(sum(abs(df[["diff"]]) <= 4*err) / nrow(df) > 0.98)

  # disagreements with WienR
  df_bad <- df[abs(df[["diff"]]) > 4*err, ]

  # small-time and large-time in fddm are equal
  fddm_sm <- da_dfddm(rt = df_bad[["rt"]], response = "lower",
                      a = df_bad[["a"]], v = df_bad[["v"]], w = df_bad[["w"]],
                      sv = df_bad[["sv"]], t0 = 0, sl_thresh = Inf,
                      err_tol = err)
  fddm_lg <- da_dfddm(rt = df_bad[["rt"]], response = "lower",
                      a = df_bad[["a"]], v = df_bad[["v"]], w = df_bad[["w"]],
                      sv = df_bad[["sv"]], t0 = 0, sl_thresh = -1,
                      err_tol = err)
  expect_true(sum(abs(fddm_sm - fddm_lg) <= 2*err) / nrow(df_bad) == 1)

  # disagreements are small relative to the calculation
  # WienR precision not guaranteed for sv > 0
  df_bad <- df_bad[df_bad[["sv"]] == 0, ]
  expect_true(all(abs(df_bad[["diff"]]) / abs(df_bad[["res_WienR"]]) <=
                  err*4e-2))



  ### dv ###
  df[["res_WienR"]] <- dvWienerPDF(t = rt, response = resp, a = a, v = v, w = w,
                                  t0 = t0, sv = sv)[["deriv"]]
  df[["res_fddm"]] <- dv_dfddm(rt = rt, response = resp, a = a, v = v, t0 = t0,
                              w = w, sv = sv, err_tol = err)
  df[["diff"]] <- df[["res_WienR"]] - df[["res_fddm"]]

  # agreement with WienR
  expect_true(sum(abs(df[["diff"]]) <= 4*err) / nrow(df) > 0.99)

  # disagreements with WienR
  df_bad <- df[abs(df[["diff"]]) > 4*err, ]

  # disagreements only for sv > 0 (WienR precision not guaranteed for sv > 0)
  expect_true(all(df_bad[["sv"]] > 0))

  # small-time and large-time in fddm are equal
  fddm_sm <- dv_dfddm(rt = df_bad[["rt"]], response = "lower",
                      a = df_bad[["a"]], v = df_bad[["v"]], w = df_bad[["w"]],
                      sv = df_bad[["sv"]], t0 = 0, sl_thresh = Inf,
                      err_tol = err)
  fddm_lg <- dv_dfddm(rt = df_bad[["rt"]], response = "lower",
                      a = df_bad[["a"]], v = df_bad[["v"]], w = df_bad[["w"]],
                      sv = df_bad[["sv"]], t0 = 0, sl_thresh = -1,
                      err_tol = err)
  expect_true(sum(abs(fddm_sm - fddm_lg) <= 2*err) / nrow(df_bad) == 1)




  ### dw ###
  df[["res_WienR"]] <- dwWienerPDF(t = rt, response = resp, a = a, v = v, w = w,
                                  t0 = t0, sv = sv)[["deriv"]]
  df[["res_fddm"]] <- dw_dfddm(rt = rt, response = resp, a = a, v = v, t0 = t0,
                              w = w, sv = sv, err_tol = err)
  df[["diff"]] <- df[["res_WienR"]] - df[["res_fddm"]]

  # agreement with WienR
  expect_true(sum(abs(df[["diff"]]) <= 4*err) / nrow(df) > 0.98)

  # disagreements with WienR
  df_bad <- df[abs(df[["diff"]]) > 4*err, ]

  # small-time and large-time in fddm are equal
  fddm_sm <- dw_dfddm(rt = df_bad[["rt"]], response = "lower",
                      a = df_bad[["a"]], v = df_bad[["v"]], w = df_bad[["w"]],
                      sv = df_bad[["sv"]], t0 = 0, sl_thresh = -1,
                      err_tol = err)
  fddm_lg <- dw_dfddm(rt = df_bad[["rt"]], response = "lower",
                      a = df_bad[["a"]], v = df_bad[["v"]], w = df_bad[["w"]],
                      sv = df_bad[["sv"]], t0 = 0, sl_thresh = Inf,
                      err_tol = err)
  expect_true(sum(abs(fddm_sm - fddm_lg) <= 2*err) / nrow(df_bad) == 1)



  ### dsv ###
  # (redefine parameter space so that sv > 0)
  if (identical(Sys.getenv("NOT_CRAN"), "true")) { # not on CRAN
    # These take a while to run
    # RT <- c(0.001, 0.01, seq(0.1, 10, by = 0.1), seq(15, 30, by = 5))
    # A <- c(0.25, seq(0.5, 5, by = 0.5))
    RT <- c(0.001, 0.1, 1, 2, 3, 4, 5, 10, 30)
    A <- c(0.25, 0.5, 1, 2.5, 5)
    V <- c(-5, -2, 0, 2, 5)
    W <- c(0.2, 0.5, 0.8)
    SV <- c(0.5, 1, 1.5)
  } else { # on CRAN
    RT <- c(0.001, 0.1, 1, 10)
    A <- c(0.5, 1, 5)
    V <- c(-5, 0, 5)
    W <- c(0.2, 0.5, 0.8)
    SV <- c(0.5, 1.5)
  }
  t0 <- 0
  err <- 1e-12 # default setting from WienR

  nRT <- length(RT)
  nA <- length(A)
  nV <- length(V)
  nW <- length(W)
  nSV <- length(SV)
  N <- nRT * nA * nV * nW * nSV

  resp <- rep("lower", N)
  rt <- rep(RT, each = nSV * nW * nV * nA, times = 1)
  a <-  rep(A,  each = nSV * nW * nV, times = nRT)
  v <-  rep(V,  each = nSV * nW, times = nRT * nA)
  w <-  rep(W,  each = nSV, times = nRT * nA * nV)
  sv <- rep(SV, each = 1, times = nRT * nA * nV * nW)

  df <- data.frame(
    rt = rt - t0,
    a = a,
    v = v,
    w = w,
    sv = sv,
    res_WienR = dsvWienerPDF(t = rt, response = resp, a = a, v = v, w = w,
                            t0 = t0, sv = sv)[["deriv"]],
    res_fddm = dsv_dfddm(rt = rt, response = resp, a = a, v = v, t0 = t0, w = w,
                        sv = sv, err_tol = err),
    diff = numeric(N)
  )
  df[["diff"]] <- df[["res_WienR"]] - df[["res_fddm"]]

  # agreement with WienR
  expect_true(sum(abs(df[["diff"]]) <= 4*err) / nrow(df) > 0.99)

  # disagreements with WienR
  df_bad <- df[abs(df[["diff"]]) > 4*err, ]

  # disagreements only for sv > 0 (WienR precision not guaranteed for sv > 0)
  expect_true(all(df_bad[["sv"]] > 0))

  # small-time and large-time in fddm are equal
  fddm_sm <- dsv_dfddm(rt = df_bad[["rt"]], response = "lower",
                      a = df_bad[["a"]], v = df_bad[["v"]], w = df_bad[["w"]],
                      sv = df_bad[["sv"]], t0 = 0, sl_thresh = Inf,
                      err_tol = err)
  fddm_lg <- dsv_dfddm(rt = df_bad[["rt"]], response = "lower",
                      a = df_bad[["a"]], v = df_bad[["v"]], w = df_bad[["w"]],
                      sv = df_bad[["sv"]], t0 = 0, sl_thresh = -1,
                      err_tol = err)
  expect_true(sum(abs(fddm_sm - fddm_lg) <= 2*err) / nrow(df_bad) == 1)
})
