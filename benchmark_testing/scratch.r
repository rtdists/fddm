library("reshape2")
library("ggplot2")


kl_Nav <- function(t, eps) {
  if (pi * t * eps < 1) {
    kl <- sqrt(-2 * log(pi*t*eps) / (pi*pi*t))
    kl <- ceiling(max(kl, 1 / (pi*sqrt(t))))
  } else {
    kl <- ceiling(1 / (pi * sqrt(t)))
  }
  return(kl)
}

lt <- seq(0.01, 5, by = 0.01)
leps <- 1e-21*10^(1:20)
kl <- data.frame(matrix(nrow = length(lt), ncol = length(leps)+1))
colnames(kl) <- c("t", leps)
kl[,1] <- lt
for (ti in 1:length(lt)) {
  for (epsi in 1:length(leps)) {
    kl[ti, epsi+1] <- kl_Nav(lt[ti], leps[epsi])
  }
}
kl_melt <- melt(kl, measure.vars = -1, variable.name = "eps", value.name = "kl")

ggplot(kl_melt, aes(x = t, y = eps)) +
  geom_raster(aes(fill = kl)) +
  scale_fill_continuous(limits=c(1, 5), breaks=1:5) +
  labs(title = "Number of terms required",
       subtitle = "Large-time summation, Navarro 2009",
       x = "Scaled Time: t/a^2 (s)", y = "Error Tolerance") +
  theme_bw() +
  theme(plot.title = element_text(size = 23),
        plot.subtitle = element_text(size = 16),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 16))




ks_Kes <- function(t, w, eps) {
  u_eps = min(-1.0, log(2 * pi * t*t * eps*eps))
  arg = -t * (u_eps - sqrt(-2 * u_eps - 2))
  k1 = (sqrt(2 * t) - w)/2
  if (arg > 0) {
    k2 = (sqrt(arg) - w) / 2;
    k = ceiling(max(k1, k2));
  }
  else {
    k = ceiling(k1);
  }
  return(k)
}

st <- seq(0.1, 5, by = 0.1)
seps <- 1e-21*10^(1:20)
sw <- seq(0.1, 0.9, by = 0.1)
ks <- data.frame(matrix(nrow = length(st)*length(seps)*length(sw), ncol = 4))
colnames(ks) <- c("t", "w", "eps", "ks")
idx <- 1
for (ti in 1:length(st)) {
  for (wi in 1:length(sw)) {
    for (epsi in 1:length(seps)) {
      ks[idx, 1] <- st[ti]
      ks[idx, 2] <- sw[wi]
      ks[idx, 3] <- seps[epsi]
      ks[idx, 4] <- ks_Kes(st[ti], sw[wi], seps[epsi])
      idx = idx + 1
    }
  }
}

ks
ks5 <- subset(ks, w == 0.5)[,c(1,3,4)]
rownames(ks5) <- NULL
melt(ks, measure.vars = 4, value.name = "ks")


st <- seq(0.01, 5, by = 0.01)
seps <- 1e-21*10^(1:20)
ks <- data.frame(matrix(nrow = length(st), ncol = length(seps)+1))
colnames(ks) <- c("t", leps)
ks[,1] <- st
for (ti in 1:length(st)) {
  for (epsi in 1:length(seps)) {
    ks[ti, epsi+1] <- ks_Kes(st[ti], 0.9, seps[epsi])
  }
}
ks_melt <- melt(ks, measure.vars = -1, variable.name = "eps", value.name = "ks")

ggplot(ks_melt, aes(x = t, y = eps, fill = ks)) +
  geom_raster() +
  scale_fill_continuous(limits=c(1, 5), breaks=1:5) +
  labs(title = "Number of terms required",
       subtitle = "Small-time summation, Kesselmeier 2014",
       x = "Scaled Time: t/a^2 (s)", y = "Error Tolerance") +
  theme_bw() +
  theme(plot.title = element_text(size = 23),
        plot.subtitle = element_text(size = 16),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 16))
