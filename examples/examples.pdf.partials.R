# derivative with respect to rt (response time)
dt_dfddm(rt = 1.2, response = "lower", a = 1, v = -1, t0 = 0.3, w = 0.5,
         sv = 0.4, err_tol = 1e-6)

# derivative with respect to t0 (non-decision time)
dt0_dfddm(rt = 1.2, response = "lower", a = 1, v = -1, t0 = 0.3, w = 0.5,
          sv = 0.4, err_tol = 1e-6)

# derivative with respect to a (threshold separation)
da_dfddm(rt = 1.2, response = "lower", a = 1, v = -1, t0 = 0.3, w = 0.5,
         sv = 0.4, err_tol = 1e-6)

# derivative with respect to v (drift rate)
dv_dfddm(rt = 1.2, response = "lower", a = 1, v = -1, t0 = 0.3, w = 0.5,
         sv = 0.4, err_tol = 1e-6)

# derivative with respect to w (relative starting point)
dw_dfddm(rt = 1.2, response = "lower", a = 1, v = -1, t0 = 0.3, w = 0.5,
         sv = 0.4, err_tol = 1e-6)

# derivative with respect to sv (inter-trial variability in the drift rate)
dsv_dfddm(rt = 1.2, response = "lower", a = 1, v = -1, t0 = 0.3, w = 0.5,
          sv = 0.4, err_tol = 1e-6)
