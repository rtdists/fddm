library("devtools")
load_all(recompile = TRUE)
library("fddm")

t <- c(0.001, 0.005, 0.01, 0.05, 0.1,
       seq(0.5, 3, by = 0.5))
fddm_fast <- dfddm_fast(rt = t, a = 1, v = -.4, t0 = 0,
                       w = 0.5, err_tol = 1e-6)
fs_Fos_17 <- dfddm(rt = t, response = 0, a = 1,
                  v = -0.4, t0 = 0, w = 0.5,
                  log = FALSE, n_terms_small = "Foster",
                  summation_small = "2017", scale = "small",
                  err_tol = 1e-5)
max(abs(fddm_fast - fs_Fos_17)) < 2*1e-5

ks_Kes_t <- function(t, w, eps) {
  k <- rep(0, length(t))
  for (i in 1:length(t)) {
    u_eps = min(-1.0, log(2 * pi * t[i]*t[i] * eps*eps))
    arg = -t[i] * (u_eps - sqrt(-2 * u_eps - 2))
    k1 = (sqrt(2 * t[i]) - w)/2
    if (arg > 0) {
      k2 = (sqrt(arg) - w) / 2;
      k[i] = ceiling(max(k1, k2));
    }
    else {
      k[i] = ceiling(k1);
    }
  }
  return(k)
}


kl_Nav_tif <- function(t, eps) {
  kl <- rep(0, length(t))
  for (i in 1:length(t)) {
    kl[i] <- ceiling(sqrt(-2 * log(pi*t[i]*eps) / (pi*pi*t[i])))
  }
  return(kl)
}
kl_Nav_tel <- function(t, eps) {
  kl <- rep(0, length(t))
  for (i in 1:length(t)) {
    kl[i] <- ceiling(1 / (pi * sqrt(t[i])))
  }
  return(kl)
}
kl_Nav_t <- function(t, eps) {
  kl <- rep(0, length(t))
  for (i in 1:length(t)) {
    if (eps < (1 / (pi*t[i]))) { # pi * eps * t[i] < 1
      kl[i] <- sqrt(-2 * log(pi*t[i]*eps) / (pi*pi*t[i]))
      kl[i] <- ceiling(max(kl[i], 1 / (pi*sqrt(t[i]))))
    } else {
      print("else")
      kl[i] <- ceiling(1 / (pi * sqrt(t[i])))
    }
  }
  return(kl)
}


t <- c(0.001, 0.005, 0.01, 0.05, 0.1,
       seq(0.5, 3, by = 0.5))
       # seq(4, 10, by = 1),
       # seq(12.5, 20, by = 2.5),
       # seq(25, 30, by = 5))
eps <- 1e-6

plot(t, kl_Nav_tif(t, eps), type = "l", col = "blue", ylim = c(0, 10))
lines(t, kl_Nav_tel(t, eps), type = "l", col = "red")
points(t, kl_Nav_t(t, eps), pch = 4, col = "purple")
lines(t, ks_Kes_t(t, 0.5, eps), type = "l", lty = "dotted", col = "darkgreen")
