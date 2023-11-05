# prepare data
data(med_dec, package = "fddm")
med_dec <- med_dec[which(med_dec[["rt"]] >= 0), ]
p1 <- med_dec[med_dec[["id"]] == 2 & med_dec[["group"]] == "experienced", ]

fit0 <- ddm(rt + response ~ 0 + classification:difficulty,
            boundary = 2, ndt = 0.39, data = p1)

fit1 <- ddm(rt + response ~ 0 + classification:difficulty,
            data = p1)

fit2 <- ddm(rt + response ~ 0 + classification:difficulty,
            bias = ~ 1, data = p1)

fit3 <- ddm(rt + response ~ 0 + classification:difficulty,
            boundary = 2, bias = ~ 1, data = p1)

fit4 <- ddm(rt + response ~ 0 + classification:difficulty,
            ndt = 0.39, data = p1)

fit5 <- ddm(rt + response ~ 0 + classification:difficulty,
            ndt = 0.39, bias = ~1, data = p1)

fit6 <- ddm(rt + response ~ 0 + classification:difficulty,
            boundary = ~ difficulty, ndt = ~ difficulty,
            bias = ~ 1, sv = ~1, data = p1,
            args_optim = list(control = list(eval.max = 600,
                                             iter.max = 600)))

test_that("ddm objects exist", {
  expect_s3_class(fit0, "ddm")
  expect_s3_class(fit1, "ddm")
  expect_s3_class(fit2, "ddm")
  expect_s3_class(fit3, "ddm")
  expect_s3_class(fit4, "ddm")
  expect_s3_class(fit5, "ddm")
  expect_s3_class(fit6, "ddm")
})

test_that("printing shows coefficients", {
  expect_output(print(fit0),
                regexp = "Fixed: boundary = 2, ndt = 0.39, bias = 0.5, sv = 0",
                fixed = TRUE)
  expect_output(print(fit0),
                regexp = "drift coefficients (identity link):", fixed = TRUE)

  expect_output(print(fit1),
                regexp = "DDM fit with 3 estimated and 2 fixed distributional parameters.\nFixed: bias = 0.5, sv = 0",
                fixed = TRUE)
  expect_output(print(fit1),
                regexp = "boundary coefficients (identity link):", fixed = TRUE)
  expect_output(print(fit1),
                regexp = "ndt coefficients (identity link):", fixed = TRUE)

  expect_output(print(fit2),
                regexp = "DDM fit with 4 estimated and 1 fixed distributional parameters.\nFixed: sv = 0",
                fixed = TRUE)
  expect_output(print(fit2),
                regexp = "bias coefficients (identity link):", fixed = TRUE)
  expect_output(print(fit2),
                regexp = "ndt coefficients (identity link):", fixed = TRUE)

  expect_output(print(fit6),
                regexp = "DDM fit with 5 estimated and 0 fixed distributional parameters.",
                fixed = TRUE)
  expect_output(print(fit6),
                regexp = "ndt coefficients (identity link):", fixed = TRUE)
  expect_output(print(fit6),
                regexp = "bias coefficients (identity link):", fixed = TRUE)
  expect_output(print(fit6),
                regexp = "sv coefficients (identity link):", fixed = TRUE)
})

test_that("correct distributional parameters are estimated", {
  expect_equal(fit0$dpar, "drift")
  expect_equal(fit1$dpar, c("drift", "boundary", "ndt"))
  expect_equal(fit2$dpar, c("drift", "boundary", "ndt", "bias"))
  expect_equal(fit3$dpar, c("drift", "ndt", "bias"))
  expect_equal(fit4$dpar, c("drift", "boundary"))
  expect_equal(fit5$dpar, c("drift", "boundary", "bias"))
  expect_equal(fit6$dpar, c("drift", "boundary", "ndt", "bias", "sv"))
})


test_that("rank deficient formulas work", {
  expect_warning(
    fit_rd <- ddm(rt + response ~ classification:difficulty, data = p1),
    regexp = "rank deficient")
  expect_equal(logLik(fit1), logLik(fit_rd))
})

# check for user supplied initial values and bounds (all combos)
test_that("initial values and bounds are checked properly", {
  # incorrect length of initial values/bounds
  expect_error(ddm(rt + response ~ 0 + classification:difficulty, data = p1,
                   args_optim = list(init = c(1, 1))))
  expect_error(ddm(rt + response ~ 0 + classification:difficulty, data = p1,
                   args_optim = list(lo_bds = c(0, 0),
                                     up_bds = c(2, 2))))
  # user supplied initial value outside bounds
  expect_warning(ddm(rt + response ~ 0 + classification:difficulty, data = p1,
                     args_optim = list(init = c(3, 3, 3, 3, 3, 3),
                                       lo_bds = c(0, 0, 0, 0, 0, 0),
                                       up_bds = c(2, 2, 2, 2, 2, 2))))
  # user supplied initial value outside generated bounds
  expect_warning(ddm(rt + response ~ 0 + classification:difficulty,
                     data = p1,
                     args_optim = list(init = c(1, 1, 1, 1, -1, -1))))
  # user supplied bounds do not include generated initial value
  expect_warning(ddm(rt + response ~ 0 + classification:difficulty, data = p1,
                     args_optim = list(lo_bds = c(50, 50, 50, 50, 50, 50),
                                       up_bds = c(70, 70, 70, 70, 50, 50))))
})
