library("Rcpp")
library("microbenchmark")

Rcpp::sourceCpp(file = "scratch.cpp")
t <- 1
a <- 1
v <- -1
w <- .5
sv <- .1
k <- 100





t <- seq(0, 3, by = 0.01)
ks_nav <- 1
eps_l_nav <- function(t, ks_nav) {
  return (1/(pi*t) * exp(-pi*pi * t * ks_nav*ks_nav / 2))
}

eps_s_nav <- function(t, ks_nav) {
  return (1/(2 * sqrt(2*pi*t)) * exp(-(ks_nav-2)*(ks_nav-2)/(2*t)))
}

ks

plot(t, eps_l_nav(t, 1), type = "l", col = "blue", ylim = c(0, 1e-6))
  lines(t, eps_l_nav(t, 2), type = "l", col = "red")
  lines(t, eps_l_nav(t, 3), type = "l", col = "darkgreen")

plot(t, eps_s_nav(t, 5), type = "l", col = "blue", ylim = c(0, 1e-6))
plot(t, eps_s_nav(t, 2), type = "l", col = "red")
plot(t, eps_s_nav(t, 3), type = "l", col = "darkgreen")


kl_nav <- function(t, eps) {
  sqrt(-2*log(pi*t*eps) / (pi*pi * t))
}

t <- 1
a <- 5
v <- 5
w <- 0.8
sv <- 0.0
mult <- exp((sv*sv * a*a * w*w - 2 * v * a * w - v*v * t) / (2 + 2 * sv*sv * t)) / (a*a * sqrt(1 + sv*sv * t))
kl_nav(t/(a*a), 1e-6/mult)

ks_nav <- function(t, eps) {
  2 + sqrt(-2*t*log(2*eps*sqrt(2*pi*t)))
}

t <- seq(0, 3, by = 0.01)
eps <- 1e-6

plot(t, kl_nav(t, eps), type = "l", col = "blue")
  lines(t, ks_nav(t, eps), type = "l", col = "red")







t <- 100
a <- 1
w <- 0.5
eps <- 1e-8

ks_eps(t, a, w, eps)
2*ks_Gon(t, w, eps)+1
ks_Nav(t, w, eps)


t <- 10^(seq(-4, 2, by = 0.5))
eps <- 10^(-seq_len(12))
ks <- matrix(nrow = length(eps), ncol = length(t))
for (i in 1:length(eps)) {
  for (j in 1:length(t)) {
    ks[i, j] <- ks_Nav(t[j], w, eps[i]) - (2*ks_Gon(t[j], w, eps[i]) + 1)
  }
}


ks

library("pheatmap")

pheatmap(ks, display_numbers = TRUE)

rep(c(2, 5, 7, 9, 12), 2)
rep(c(1, 2), each = 5)
factor(rep(1:5, each = 2))






mbm <- microbenchmark(
  ks_1 = ks_eps_17(1, 1, 0.5, 1e-3),
  ks_2 = ks_eps_172(1, 1, 0.5, 1e-3),
  times = 100000
)
summary(mbm)



t <- 10^(seq(-4, 2, by = 0.05))
a <- 1
eps <- 10^(-seq_len(12))
w <- seq(0, 1, by = 0.1)

FuncNames <- c("t", "eps", "w",
               "kb_S4", "kb_S7", "kb_G", "kb_N",
               "ks_S4", "ks_S7", "ks_G", "ks_N", "kl_N")
kb <- data.frame(matrix(nrow = length(t)*length(eps)*length(w),
                        ncol = length(FuncNames)))
colnames(kb) <- FuncNames

for (i in 1:length(t)) {
  for (j in 1:length(eps)) {
    for (k in 1:length(w)) {
      ijk <- (i-1)*length(eps)*length(w) + (j-1)*length(w) + k
      kb$t[ijk] <- i
      kb$eps[ijk] <- j
      kb$w[ijk] <- k
      kb$kb_S4[ijk] <- min(ks_eps_14(t[i], a, w[k], eps[j]), kl_Nav(t[i]/a/a, w[k], eps[j]))
      kb$kb_S7[ijk] <- min(ks_eps_17(t[i], a, w[k], eps[j]), kl_Nav(t[i]/a/a, w[k], eps[j]))
      kb$kb_G[ijk] <- min(ks_Gon(t[i]/a/a, w[k], eps[j]), kl_Nav(t[i]/a/a, w[k], eps[j]))
      kb$kb_N[ijk] <- min(ks_Nav(t[i]/a/a, w[k], eps[j]), kl_Nav(t[i]/a/a, w[k], eps[j]))
      kb$ks_S4[ijk] <- ks_eps_14(t[i], a, w[k], eps[j])
      kb$ks_S7[ijk] <- ks_eps_17(t[i], a, w[k], eps[j])
      kb$ks_G[ijk] <- ks_Gon(t[i]/a/a, w[k], eps[j])
      kb$ks_N[ijk] <- ks_Nav(t[i]/a/a, w[k], eps[j])
      kb$kl_N[ijk] <- kl_Nav(t[i]/a/a, w[k], eps[j])
    }
  }
}



kb$ks_S4 == kb$ks_S7
summary(kb$kb_S4)
summary(kb$kb_S7)
summary(kb$kb_G)
summary(kb$kb_N)



length(which(kb[["ks_S7"]] <  kb[["ks_S4"]]))
length(which(kb[["ks_S7"]] == kb[["ks_S4"]]))
length(which(kb[["ks_S7"]] >  kb[["ks_S4"]]))
t[81]
eps[6]
w[1]
max(kb[kb[["ks_S7"]] >  kb[["ks_G"]], ]$eps)


kN <- 6
length(-floor((kN-1)/2):ceiling((kN-1)/2))


library("ggplot2")
ggplot(kb) +
  geom_tile(aes(x = t, y = eps, fill = w))

ggplot(subset(kb, w == 6)) +
  geom_tile(aes(x = t, y = eps, fill = kb_S7))




plot(t, ks_e, type = "l", lwd = 4, col = "green", log = "x")
  lines(t, ks_G, type = "l", col = "blue")
  lines(t, ks_N, type = "l", col = "red")
  lines(t, kl_N, type = "l", lwd = 2, col = "brown")

plot(t, pmin(ks_e, kl_N), type = "l", lwd = 4, col = "green", log = "x", ylim = c(1, 5))
  lines(t, pmin(ks_G, kl_N), type = "l", col = "blue")
  lines(t, pmin(ks_N, kl_N), type = "l", col = "red")







t <- 81
a <- 1
w <- 1
gamma <- -a*a / (2 * t)
j <- 1:10

sqt <- sqrt(t)
jp <- 0.5 * (-w + sqt/a)
jn <- 0.5 * (-w - sqt/a)
jp
jn
abs(jp - jn)
abs(abs(jp) - abs(jn))

pterm <- function(j, gamma, w) {
  return((w + 2 * j) * exp(gamma * (w + 2 * j) * (w + 2 * j)))
}
nterm <- function(j, gamma, w) {
  return(abs((w - 2 * j) * exp(gamma * (w - 2 * j) * (w - 2 * j))))
}
jvec <- seq(j[1], j[length(j)], by = 0.1)
jpvec <- pterm(jvec, gamma, w)
jnvec <- nterm(jvec, gamma, w)
plot(jvec, jpvec, type = "l", col = "blue",
  xlim = c(0, j[length(j)]),
  ylim = c(min(jpvec, jnvec), max(jpvec, jnvec)))
  lines(jvec, jnvec, col = "red")
  points(j, pterm(j, gamma, w), col = "blue", pch = 1, cex = 2)
  points(j, nterm(j, gamma, w), col = "red", pch = 16, cex = 1.5)
  abline(v = abs(jp), lty = "dashed", col = "blue")
  abline(v = abs(jn), lty = "dashed", col = "red")

jpn <- mean(c(jp, abs(jn)))

pterm(jpn, gamma, w) - nterm(jpn, gamma, w)

jvec <- seq(j[1], j[length(j)], by = 0.1)
jpvec <- pterm(jvec, gamma, w)
jnvec <- nterm(jvec, gamma, w)
jvec[which.min(abs(jpvec - jnvec))] == mean(c(abs(0.5 * (-w + sqt/a)), abs(0.5 * (-w - sqt/a))))
jvec[which.min(abs(jpvec - jnvec))] == sqrt(t)/(2*a)

pterm(5, gamma, w) == nterm(6, gamma, w)
