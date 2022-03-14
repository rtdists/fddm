# minimal example
pfddm(rt = 1.2, response = "lower", a = 1, v = -1, t0 = 0.3)

# example with all function parameters set to default or a practical value
pfddm(rt = c(1.2, 0.9, 1.1, 1.4, 0.8, 1.3),
      response = c("lower", "upper", "upper", "lower", "upper", "lower"),
      a = 1, v = -1, t0 = 0.2, w = 0.5, sv = 0, sigma = 1,
      err_tol = 1e-6, log = FALSE, method = "Mills")

# example of mismatched input lengths
pfddm(rt = c(1.2, 0.9, 1.1, 1.4, 0.8, 1.3),
      response = c("lower", "upper", "upper", "lower", "upper", "lower"),
      a = c(1, 3), v = c(-2, 2, 2, -2, 2, -2),
      t0 = 0.3, w = c(0.4, 0.5, 0.6), sv = 0.9,
      err_tol = 1e-10, log = FALSE, method = "NCDF")

# example with Wiener diffusion coefficient (sigma) not equal to 1
pfddm(rt = c(1.2, 0.9, 1.1, 1.4, 0.8, 1.3),
      response = c("lower", "upper", "upper", "lower", "upper", "lower"),
      a = 1, v = -1, t0 = 0.3, w = 0.5, sv = 0, sigma = 0.1,
      err_tol = 1e-10, log = TRUE, method = "Mills")


### examples of different response inputs

# integer
resp_int <- as.integer(c(1, 2, 2, 1, 2, 1))
pfddm(rt = c(1.2, 0.9, 1.1, 1.4, 0.8, 1.3), response = resp_int,
      a = 1, v = -2, t0 = 0.3, w = 0.5, sv = 0.1,
      err_tol = 1e-10, log = FALSE, method = "Mills")

# double
resp_dbl <- as.double(c(1, 2, 2, 1, 2, 1))
pfddm(rt = c(1.2, 0.9, 1.1, 1.4, 0.8, 1.3), response = resp_dbl,
      a = 1, v = -2, t0 = 0.3, w = 0.5, sv = 0.1,
      err_tol = 1e-10, log = FALSE, method = "Mills")

# factor (first level is mapped to "lower")
days <- c("Monday", "Friday", "Friday", "Monday", "Friday", "Monday")
resp_fac <- factor(days, levels = c("Monday", "Friday"))
pfddm(rt = c(1.2, 0.9, 1.1, 1.4, 0.8, 1.3), response = resp_fac,
      a = 1, v = -2, t0 = 0.3, w = 0.5, sv = 0.1,
      err_tol = 1e-10, log = FALSE, method = "Mills")

# string
resp_str <- c("lower", "upper", "upper", "lower", "upper", "lower")
pfddm(rt = c(1.2, 0.9, 1.1, 1.4, 0.8, 1.3), response = resp_str,
      a = 1, v = -2, t0 = 0.3, w = 0.5, sv = 0.1,
      err_tol = 1e-10, log = FALSE, method = "Mills")

# logical
resp_log <- c(FALSE, TRUE, TRUE, FALSE, TRUE, FALSE)
pfddm(rt = c(1.2, 0.9, 1.1, 1.4, 0.8, 1.3), response = resp_log,
      a = 1, v = -2, t0 = 0.3, w = 0.5, sv = 0.1,
      err_tol = 1e-10, log = FALSE, method = "Mills")
