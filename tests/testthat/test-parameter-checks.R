context("Testing the parameter checks for consistency")
# this pertains to dfddm() and pfddm() as they share the functionality

### Input checking
test_that("Input checking", {
  
  # rt (dfddm and pfddm have different behavior for t -> +Inf)
  expect_equal(
    dfddm(rt = c(0, -1, Inf, -Inf), response = 1, a = 1, v = -1, t0 = 0,
          w = 0.5, sv = 0, sigma = 1, log = 0, err_tol = 1e-6),
    c(0, 0, 0, 0)
  )
  expect_equal(
    dfddm(rt = c(0, -1, Inf, -Inf), response = 1, a = 1, v = -1, t0 = 0,
          w = 0.5, sv = 0, sigma = 1, log = 1, err_tol = 1e-6),
    c(-Inf, -Inf, -Inf, -Inf)
  )
  expect_equal(
    pfddm(rt = c(-Inf, Inf), response = 1, a = 1, v = -1, t0 = 0,
          w = 0.5, sv = 0, sigma = 1, log = 0, err_tol = 1e-6),
    c(0, 0.7310586),
    tolerance = 1e-6
  )
  expect_equal(
    pfddm(rt = c(-Inf, Inf), response = 1, a = 1, v = -1, t0 = 0,
          w = 0.5, sv = 0, sigma = 1, log = 1, err_tol = 1e-6),
    c(-Inf, log(0.7310586)),
    tolerance = 1e-6
  )
  expect_true(all(c(
    is.nan(dfddm(rt = NaN, response = 1, a = 1, v = -1, t0 = 0, w = 0.5,
                 sv = 0, sigma = 1, log = 0, err_tol = 1e-6)),
    is.na(dfddm(rt = NA, response = 1, a = 1, v = -1, t0 = 0, w = 0.5,
                sv = 0, sigma = 1, log = 0, err_tol = 1e-6))
  )))
  expect_true(all(c(
    is.nan(dfddm(rt = NaN, response = 1, a = 1, v = -1, t0 = 0, w = 0.5,
                 sv = 0, sigma = 1, log = 1, err_tol = 1e-6)),
    is.na(dfddm(rt = NA, response = 1, a = 1, v = -1, t0 = 0, w = 0.5,
                sv = 0, sigma = 1, log = 1, err_tol = 1e-6))
  )))
  expect_warning(expect_equal(
    length(dfddm(rt = numeric(), response = 1, a = 1, v = -1, t0 = 0, w = 0.5,
                 sv = 0, sigma = 1, log = 0, err_tol = 1e-6)),
    0
  ))
  
  # response
  expect_warning(expect_equal(
    dfddm(rt = 1, response = c(3, -1, Inf, -Inf), a = 1, v = -1, t0 = 0,
          w = 0.5, sv = 0, sigma = 1, log = 0, err_tol = 1e-6),
    c(0, 0, 0, 0)
  ))
  expect_warning(expect_true(all(c(
    is.nan(dfddm(rt = 1, response = NaN, a = 1, v = -1, t0 = 0, w = 0.5,
                 sv = 0, sigma = 1, log = 0, err_tol = 1e-6)),
    is.na(dfddm(rt = 1, response = NA, a = 1, v = -1, t0 = 0, w = 0.5,
                sv = 0, sigma = 1, log = 0, err_tol = 1e-6))
  ))))
  expect_warning(expect_equal(
    dfddm(rt = 1, response = c("l", "u", Inf, -Inf, NaN, NA), a = 1, v = -1,
          t0 = 0, w = 0.5, sv = 0, sigma = 1, log = 0, err_tol = 1e-6)[1:4],
    c(dfddm(rt = 1, response = c("l", "u"), a = 1, v = -1, t0 = 0, w = 0.5,
            sv = 0, sigma = 1, log = 0, err_tol = 1e-6), 0, 0)
  ))
  expect_warning(expect_true(all(c(
    is.nan(dfddm(rt = 1, response = c("l", "u", Inf, -Inf, NaN, NA), a = 1,
                 v = -1, t0 = 0, w = 0.5, sv = 0, sigma = 1, log = 0,
                 err_tol = 1e-6)[5]),
    is.na(dfddm(rt = 1, response = c("l", "u", Inf, -Inf, NaN, NA), a = 1,
                v = -1, t0 = 0, w = 0.5, sv = 0, sigma = 1, log = 0,
                err_tol = 1e-6)[6])
  ))))
  expect_warning(expect_equal(
    length(dfddm(rt = 1, response = numeric(), a = 1, v = -1, t0 = 0, w = 0.5,
                 sv = 0, sigma = 1, log = 0, err_tol = 1e-6)),
    0
  ))
  expect_warning(expect_equal(
    length(dfddm(rt = 1, response = character(), a = 1, v = -1, t0 = 0, w = 0.5,
                 sv = 0, sigma = 1, log = 0, err_tol = 1e-6)),
    0
  ))
  expect_warning(expect_equal(
    length(dfddm(rt = 1, response = logical(), a = 1, v = -1, t0 = 0, w = 0.5,
                 sv = 0, sigma = 1, log = 0, err_tol = 1e-6)),
    0
  ))
  
  # a
  expect_true(all(c(
    is.nan(dfddm(rt = 1, response = 1, a = c(-0.4, 0, NaN, -Inf), v = -1,
                 t0 = 0, w = 0.5, sv = 0, sigma = 1, log = 0, err_tol = 1e-6)),
    is.na(dfddm(rt = 1, response = 1, a = NA, v = -10, t0 = 0,
                w = 0.5, sv = 0, sigma = 1, log = 0, err_tol = 1e-6))
  )))
  expect_equal(
    dfddm(rt = 1, response = 1, a = Inf, v = -10, t0 = 0, scale = "l",
          n_terms_small = "n", w = 0.5, sv = 0, sigma = 1, log = 0,
          err_tol = 1e-6),
    0
  )
  expect_equal(
    dfddm(rt = 1, response = 1, a = Inf, v = -10, t0 = 0, scale = "l",
          n_terms_small = "n", w = 0.5, sv = 0, sigma = 1, log = 1,
          err_tol = 1e-6),
    -Inf
  )
  expect_warning(expect_equal(
    length(dfddm(rt = 1, response = 1, a = numeric(), v = -1, t0 = 0, w = 0.5,
                 sv = 0, sigma = 1, log = 0, err_tol = 1e-6)),
    0
  ))
  
  # v
  expect_true(all(c(
    is.nan(dfddm(rt = 1, response = 1, a = 1, v = NaN, t0 = 0,
                 w = 0.5, sv = 0, sigma = 1, log = 0, err_tol = 1e-6)),
    is.na(dfddm(rt = 1, response = 1, a = 1, v = NA, t0 = 0,
                w = 0.5, sv = 0, sigma = 1, log = 0, err_tol = 1e-6))
  )))
  expect_equal(
    dfddm(rt = 1, response = 1, a = 1, v = c(Inf, -Inf), t0 = 0,
          w = 0.5, sv = 0, sigma = 1, log = 0, err_tol = 1e-6),
    c(0, 0)
  )
  expect_equal(
    dfddm(rt = 1, response = 1, a = 1, v = c(Inf, -Inf), t0 = 0,
          w = 0.5, sv = 0, sigma = 1, log = 1, err_tol = 1e-6),
    c(-Inf, -Inf)
  )
  expect_warning(expect_equal(
    length(dfddm(rt = 1, response = 1, a = 1, v = numeric(), t0 = 0, w = 0.5,
                 sv = 0, sigma = 1, log = 0, err_tol = 1e-6)),
    0
  ))
  
  # t0
  expect_true(all(c(
    is.nan(dfddm(rt = 1, response = 1, a = 1, v = -1, t0 = c(-0.25, NaN, -Inf),
                 w = 0.5, sv = 0, sigma = 1, log = 0, err_tol = 1e-6)),
    is.na(dfddm(rt = 1, response = 1, a = 1, v = -1, t0 = NA,
                w = 0.5, sv = 0, sigma = 1, log = 0, err_tol = 1e-6))
  )))
  expect_equal(
    dfddm(rt = 1, response = 1, a = 1, v = -1, t0 = Inf,
          w = 0.5, sv = 0, sigma = 1, log = 0, err_tol = 1e-6),
    0
  )
  expect_equal(
    dfddm(rt = 1, response = 1, a = 1, v = -1, t0 = Inf,
          w = 0.5, sv = 0, sigma = 1, log = 1, err_tol = 1e-6),
    -Inf
  )
  expect_warning(expect_equal(
    length(dfddm(rt = 1, response = 1, a = 1, v = -1, t0 = numeric(), w = 0.5,
                 sv = 0, sigma = 1, log = 0, err_tol = 1e-6)),
    0
  ))
  
  # w
  expect_true(all(c(
    is.nan(dfddm(rt = 1, response = 1, a = 1, v = -1, t0 = 0,
                 w = c(-0.5, 1.5, 0, 1, NaN, Inf, -Inf), sv = 0, sigma = 1,
                 log = 0, err_tol = 1e-6)),
    is.na(dfddm(rt = 1, response = 1, a = 1, v = -1, t0 = 0,
                w = NA, sv = 0, sigma = 1, log = 0, err_tol = 1e-6))
  )))
  expect_warning(expect_equal(
    length(dfddm(rt = 1, response = 1, a = 1, v = -1, t0 = 0, w = numeric(),
                 sv = 0, sigma = 1, log = 0, err_tol = 1e-6)),
    0
  ))
  
  # sv
  expect_true(all(c(
    is.nan(dfddm(rt = 1, response = 1, a = 1, v = -1, t0 = 0, w = 0.5,
                 sv = c(-1, NaN, -Inf), sigma = 1, log = 0, err_tol = 1e-6)),
    is.na(dfddm(rt = 1, response = 1, a = 1, v = -1, t0 = 0, w = 0.5,
                sv = NA, sigma = 1, log = 0, err_tol = 1e-6))
  )))
  expect_equal(
    dfddm(rt = 1, response = 1, a = 1, v = -1, t0 = 0, w = 0.5,
          sv = Inf, sigma = 1, log = 0, err_tol = 1e-6),
    0
  )
  expect_equal(
    dfddm(rt = 1, response = 1, a = 1, v = -1, t0 = 0, w = 0.5,
          sv = Inf, sigma = 1, log = 1, err_tol = 1e-6),
    -Inf
  )
  expect_warning(expect_equal(
    length(dfddm(rt = 1, response = 1, a = 1, v = -1, t0 = 0, w = 0.5,
                 sv = numeric(), sigma = 1, log = 0, err_tol = 1e-6)),
    0
  ))
  
  # sigma
  expect_true(all(c(
    is.nan(dfddm(rt = 1, response = 1, a = 1, v = -1, t0 = 0, w = 0.5, sv = 0,
                 sigma = c(-1, 0, NaN, Inf, -Inf), log = 0, err_tol = 1e-6)),
    is.na(dfddm(rt = 1, response = 1, a = 1, v = -1, t0 = 0, w = 0.5,
                sv = 0, sigma = NA, log = 0, err_tol = 1e-6))
  )))
  expect_warning(expect_equal(
    length(dfddm(rt = 1, response = 1, a = 1, v = -1, t0 = 0, w = 0.5,
                 sv = 0, sigma = numeric(), log = 0, err_tol = 1e-6)),
    0
  ))
  
  # err_tol
  expect_error(
    dfddm(rt = 1, response = 1, a = 1, v = -1, t0 = 0, w = 0.5, sv = 0,
          sigma = 1, log = 0, err_tol = c(0, -1e-6, NA, NaN, Inf, -Inf))
  )
  expect_error(
    dfddm(rt = 1, response = 1, a = 1, v = -1, t0 = 0, w = 0.5,
          sv = 0, sigma = 1, log = 0, err_tol = numeric())
  )
})
