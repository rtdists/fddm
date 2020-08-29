# minimal example
dfddm(rt = 1, response = "lower", a = 1, v = -2, t0 = 0.3)

# more practical example
dfddm(rt = c(3, 2.75, 1.4, 2, 2.3, 1.8),
      response = c("lower", "upper", "upper", "lower", "upper", "lower"),
      a = 1, v = -2, t0 = 0.3, w = 0.4, sv = 0.1,
      log = TRUE, n_terms_small = "SWSE",
      summation_small = "2017", scale = "both",
      err_tol = 1e-10)

# example of mismatched input lengths
dfddm(rt = c(3, 2.75, 1.4, 2, 2.3, 1.8),
      response = c("lower", "upper", "upper", "lower", "upper", "lower"),
      a = c(1, 3), v = c(-2, 2, 2, -2, 2, -2),
      t0 = 0.3, w = c(0.4, 0.6, 0.5), sv = 0.1,
      log = FALSE, n_terms_small = "Navarro",
      summation_small = "2017", scale = "large",
      err_tol = 1e-10)


### examples of different response inputs

# string
resp_str <- c("lower", "upper", "upper", "lower", "upper", "lower")
dfddm(rt = c(3, 2.75, 1.4, 2, 2.3, 1.8), response = resp_str,
      a = 1, v = -2, t0 = 0.3, w = 0.4, sv = 0.1,
      log = FALSE, n_terms_small = "SWSE",
      summation_small = "2017", scale = "both",
      err_tol = 1e-10)

# integer
resp_int <- as.integer(c(1, 2, 2, 1, 2, 1))
dfddm(rt = c(3, 2.75, 1.4, 2, 2.3, 1.8), response = resp_int,
      a = 1, v = -2, t0 = 0.3, w = 0.4, sv = 0.1,
      log = FALSE, n_terms_small = "SWSE",
      summation_small = "2017", scale = "both",
      err_tol = 1e-10)

# double
resp_dbl <- as.double(c(1, 2, 2, 1, 2, 1))
dfddm(rt = c(3, 2.75, 1.4, 2, 2.3, 1.8), response = resp_dbl,
      a = 1, v = -2, t0 = 0.3, w = 0.4, sv = 0.1,
      log = FALSE, n_terms_small = "SWSE",
      summation_small = "2017", scale = "both",
      err_tol = 1e-10)

# logical
resp_log <- c(FALSE, TRUE, TRUE, FALSE, TRUE, FALSE)
dfddm(rt = c(3, 2.75, 1.4, 2, 2.3, 1.8), response = resp_log,
      a = 1, v = -2, t0 = 0.3, w = 0.4, sv = 0.1,
      log = FALSE, n_terms_small = "SWSE",
      summation_small = "2017", scale = "both",
      err_tol = 1e-10)

# factor (first level is mapped to "lower")
days <- c("Monday", "Friday", "Friday", "Monday", "Friday", "Monday")
resp_fac <- factor(days, levels = c("Monday", "Friday"))
dfddm(rt = c(3, 2.75, 1.4, 2, 2.3, 1.8), response = resp_fac,
      a = 1, v = -2, t0 = 0.3, w = 0.4, sv = 0.1,
      log = FALSE, n_terms_small = "SWSE",
      summation_small = "2017", scale = "both",
      err_tol = 1e-10)
