# minimal example
dfddm(rt = 1.2, response = "lower", a = 1, v = -1, t0 = 0.3)

# example with all function parameters set to default or a practical value
dfddm(rt = c(1.2, 0.9, 1.1, 1.4, 0.8, 1.3),
      response = c("lower", "upper", "upper", "lower", "upper", "lower"),
      a = 1, v = -1, t0 = 0.2, w = 0.5, sv = 0, sigma = 1,
      err_tol = 1e-6, log = FALSE, switch_mech = "eff_rt", switch_thresh = 0.8,
      n_terms_small = "SWSE", summation_small = "2017")

# example of mismatched input lengths
dfddm(rt = c(1.2, 0.9, 1.1, 1.4, 0.8, 1.3),
      response = c("lower", "upper", "upper", "lower", "upper", "lower"),
      a = c(1, 3), v = c(-2, 2, 2, -2, 2, -2),
      t0 = 0.3, w = c(0.4, 0.5, 0.6), sv = 0.9,
      err_tol = 1e-10, log = FALSE, switch_mech = "large",
      summation_small = "2017")

# example with Wiener diffusion coefficient (sigma) not equal to 1
dfddm(rt = c(1.2, 0.9, 1.1, 1.4, 0.8, 1.3),
      response = c("lower", "upper", "upper", "lower", "upper", "lower"),
      a = 1, v = -1, t0 = 0.3, w = 0.5, sv = 0, sigma = 0.1,
      err_tol = 1e-10, log = TRUE, switch_mech = "terms_large",
      switch_thresh = 1, summation_small = "2017")


### examples of different response inputs

# integer
resp_int <- as.integer(c(1, 2, 2, 1, 2, 1))
dfddm(rt = c(1.2, 0.9, 1.1, 1.4, 0.8, 1.3), response = resp_int,
      a = 1, v = -2, t0 = 0.3, w = 0.5, sv = 0.1,
      err_tol = 1e-10, log = FALSE, switch_mech = "eff_rt",
      summation_small = "2017")

# double
resp_dbl <- as.double(c(1, 2, 2, 1, 2, 1))
dfddm(rt = c(1.2, 0.9, 1.1, 1.4, 0.8, 1.3), response = resp_dbl,
      a = 1, v = -2, t0 = 0.3, w = 0.5, sv = 0.1,
      err_tol = 1e-10, log = FALSE, switch_mech = "eff_rt",
      summation_small = "2017")

# factor (first level is mapped to "lower")
days <- c("Monday", "Friday", "Friday", "Monday", "Friday", "Monday")
resp_fac <- factor(days, levels = c("Monday", "Friday"))
dfddm(rt = c(1.2, 0.9, 1.1, 1.4, 0.8, 1.3), response = resp_fac,
      a = 1, v = -2, t0 = 0.3, w = 0.5, sv = 0.1,
      err_tol = 1e-10, log = FALSE, switch_mech = "eff_rt",
      summation_small = "2017")

# string
resp_str <- c("lower", "upper", "upper", "lower", "upper", "lower")
dfddm(rt = c(1.2, 0.9, 1.1, 1.4, 0.8, 1.3), response = resp_str,
      a = 1, v = -2, t0 = 0.3, w = 0.5, sv = 0.1,
      err_tol = 1e-10, log = FALSE, switch_mech = "eff_rt",
      summation_small = "2017")

# logical
resp_log <- c(FALSE, TRUE, TRUE, FALSE, TRUE, FALSE)
dfddm(rt = c(1.2, 0.9, 1.1, 1.4, 0.8, 1.3), response = resp_log,
      a = 1, v = -2, t0 = 0.3, w = 0.5, sv = 0.1,
      err_tol = 1e-10, log = FALSE, switch_mech = "eff_rt",
      summation_small = "2017")
