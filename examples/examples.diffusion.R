rt <- seq(0.1, 3, by = 0.1)
response <- c(rep(0, length(rt)/2), rep(1, length(rt)/2))
dfddm(rt = rt, response = response, a = 1, v = 2, t0 = 0.3, w = 0.5, sv = 0.5,
      log = FALSE, n_terms_small = "Foster", summation_small = "2017",
      scale = "small", err_tol = 0.000001)

rt <- seq(0.1, 3, by = 0.1)
response <- 0
dfddm(rt = rt, response = response, a = 1, v = 2, t0 = 0.3, w = 0.5, sv = 0.5,
      log = FALSE, n_terms_small = "Foster", summation_small = "2017",
      scale = "small", err_tol = 0.000001)

rt <- seq(0.1, 3, by = 0.1)
response <- c(rep(0, length(rt)/2), rep(1, length(rt)/2))
dfddm(rt = rt, response = response, a = 1, v = 2, t0 = 0.3, w = 0.5, sv = 0.5,
      log = TRUE, n_terms_small = "Kesselmeier", summation_small = "2014",
      scale = "both", err_tol = 0.001)

rt <- seq(0.1, 3, by = 0.1)
response <- rep(0, length(rt))
dfddm(rt = rt, response = response, a = 1.5, v = -0.4, t0 = 0.3)
