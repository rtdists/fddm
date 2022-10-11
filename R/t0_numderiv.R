# dt_wrap <- function(t, p) {
#   return(dt_dfddm(rt = t, response = "lower", v = p[1], a = p[2], t0 = 0,
#                   w = p[3], sv = p[4], err_tol = err))
# }
# numDeriv::grad(func = dt_wrap, x = c("t" = t_vec[i]), method = "Richardson",
#                # method.args = list(eps = 1.5*err, d = 0.0001, zero.tol = err,
#                #                    r = 6, v = 2, show.details = FALSE),
#                p = c(v = v, a = a, w = w, sv = sv))
