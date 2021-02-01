# This file produces the sample terms in the various timescales of the truncated
# infinite sums, as seen in Section 3 of the fddm paper:
# "Approximations to the DDM Density Functions"

library("ggplot2")
save_dir <- "paper_analysis/images/"

fl <- function(j, t, a, w) {
  return(j*sin(j*w*pi)*exp(-(j*j*pi*pi*t)/(2*a*a)))
}

fs_14 <- function(j, t, a, w) {
  return((w+2*j)*exp(-(a*a)*(w+2*j)*(w+2*j)/(2*t)))
}

fs_17 <- function(j, t, a, w) {
  gamma <- -a*a / (2*t)
  out <- rep(0, length(j))
  out[j%%2 == 1] <- -(j[j%%2 == 1]+1-w) *
    exp(gamma * (j[j%%2 == 1]+1-w)*(j[j%%2 == 1]+1-w))
  out[j%%2 == 0] <- (j[j%%2 == 0]+w) *
    exp(gamma * (j[j%%2 == 0]+w)*(j[j%%2 == 0]+w))
  return(out)
}

jl <- 0:20
js <- -10:10
tl <- 0.005
ts <- 25
a <- 1
w <- 0.69
eps <- 0.12

fl1 <- fl(jl, tl, a, w)
fs1_14 <- fs_14(js, ts, a, w)
fs1_17 <- fs_17(jl, ts, a, w)

df <- data.frame(jl = jl,
                 js = js,
                 fl = fl1,
                 fs_14 = fs1_14,
                 fs_17 = fs1_17,
                 dummy = 1,
                 dummy3 = 3)

jl2 <- seq(0, 20, by = 0.01)
js2 <- seq(-10, 10, by = 0.01)
fl2 <- fl(jl2, tl, a, w)
fs2_14 <- fs_14(js2, ts, a, w)
fs2_17 <- fs_17(jl2, ts, a, w)
df2 <- data.frame(jl = jl2,
                  js = js2,
                  fl = fl2,
                  fs_14 = fs2_14,
                  fs_17 = fs2_17,
                  eps = 0.12,
                  dummy = 2,
                  dummy3 = 3)

df3 <- data.frame(x = 13.5,
                  y = min(fl2):max(fl2),
                  dummy4 = 4)


########## value of terms plot, large-time
ggplot() +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_line(data = df2, aes(x = jl, y = fl, color = factor(dummy, levels = 1:4)), size = 1) +
  geom_point(data = df, aes(x = jl, y = fl, color = factor(dummy, levels = 1:4)), size = 4) +
  geom_line(data = df3, aes(x = x, y = y, color = factor(dummy4, levels = 1:4)), linetype = "dashed", size = 1) +
  geom_point(data = df[df[["jl"]] > 13.5,], aes(x = jl, y = fl, color = factor(dummy3, levels = 1:4)), size = 4) +
  labs(title = "Values of the terms in the large-time summation, by index",
       x = "Index of the term",
       y = "Value of the term") +
  scale_color_manual(values = c("#cb5041", "#decbc9", "blue", "blue"),
                     name = NULL,
                     breaks = 1:4,
                     labels = c("discretization", "continuous", "truncation", "post truncation")) +
  guides(color = guide_legend(override.aes = list(shape = c(16, NA, NA, 16),
                                                  linetype = c("blank", "solid", "dashed", "blank")))) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        plot.title = element_text(size = 20),
        plot.subtitle = element_text(size = 16,
                                     margin = margin(5, 5, 15, 5, "pt")),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        strip.text = element_text(size = 14),
        strip.background = element_rect(fill = "white"),
        legend.position = c(1,1),
        legend.justification = c(1,1),
        legend.box = "vertical",
        legend.direction = "vertical",
        legend.background = element_rect(fill = "transparent"),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 13),
        legend.key.width = unit(3,"line"))
ggsave(paste0(save_dir, "large_terms_eps.png"), width = 16, height = 9)



########## value of terms plot, small-time, S_14 style
ggplot() +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_line(data = df2, aes(x = js, y = fs_14, color = factor(dummy, levels = 1:3)), size = 1) +
  geom_point(data = df, aes(x = js, y = fs_14, color = factor(dummy, levels = 1:3)), size = 4) +
  geom_line(data = df2, aes(x = js, y = eps, color = factor(dummy3, levels = 1:3)), linetype = "dashed", size = 1) +
  geom_line(data = df2, aes(x = js, y = -eps, color = factor(dummy3, levels = 1:3)), linetype = "dashed", size = 1) +
  geom_point(data = df[abs(df[["fs_14"]]) < eps,], aes(x = js, y = fs_14, color = factor(dummy3, levels = 1:4)), size = 4) +
  labs(title = "Values of the terms in the small-time summation, by index",
       x = "Index of the term",
       y = "Value of the term") +
  scale_color_manual(values = c("#cb5041", "#decbc9", "blue"),
                     name = NULL,
                     breaks = c(1, 2, 3),
                     labels = c("discretization", "continuous", "error tolerance")) +
  guides(color = guide_legend(override.aes = list(shape = c(16, NA, NA),
                                                  linetype = c("blank", "solid", "dashed")))) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        plot.title = element_text(size = 20),
        plot.subtitle = element_text(size = 16,
                                     margin = margin(5, 5, 15, 5, "pt")),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        strip.text = element_text(size = 14),
        strip.background = element_rect(fill = "white"),
        legend.position = c(1,1),
        legend.justification = c(1,1),
        legend.box = "vertical",
        legend.direction = "vertical",
        legend.background = element_rect(fill = "transparent"),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 13),
        legend.key.width = unit(3,"line"))
ggsave(paste0(save_dir, "small_terms_eps_14.png"), width = 16, height = 9)


########## value of terms plot, small-time, S_17 style
ggplot() +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_point(data = df, aes(x = jl, y = fs_17, color = factor(dummy, levels = 1:3)), size = 4) +
  geom_line(data = df2, aes(x = jl, y = eps, color = factor(dummy3, levels = 1:3)), linetype = "dashed", size = 1) +
  geom_line(data = df2, aes(x = jl, y = -eps, color = factor(dummy3, levels = 1:3)), linetype = "dashed", size = 1) +
  geom_point(data = df[abs(df[["fs_17"]]) < eps,], aes(x = jl, y = fs_17, color = factor(dummy3, levels = 1:4)), size = 4) +
  labs(title = "Values of the terms in the small-time summation, by index",
       x = "Index of the term",
       y = "Value of the term") +
  scale_color_manual(values = c("#cb5041", "blue"),
                     name = NULL,
                     breaks = c(1, 3),
                     labels = c("discretization", "error tolerance")) +
  guides(color = guide_legend(override.aes = list(shape = c(16, NA),
                                                  linetype = c("blank", "dashed")))) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        plot.title = element_text(size = 20),
        plot.subtitle = element_text(size = 16,
                                     margin = margin(5, 5, 15, 5, "pt")),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        strip.text = element_text(size = 14),
        strip.background = element_rect(fill = "white"),
        legend.position = c(1,1),
        legend.justification = c(1,1),
        legend.box = "vertical",
        legend.direction = "vertical",
        legend.background = element_rect(fill = "transparent"),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 13),
        legend.key.width = unit(3,"line"))
ggsave(paste0(save_dir, "small_terms_eps_17.png"), width = 16, height = 9)
