# Figures 1, 2, 3

# This file produces the sample terms in the various timescales of the truncated
# infinite sums, as seen in Section 3 of the fddm paper:
# "Approximations to the DDM Density Functions"

library("ggplot2")
img_dir <- "paper_analysis/images/sec_3_approximations/"

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


tl <- 0.005
ts <- 25
a <- 1
w <- 0.69
eps <- 0.12

jl_c <- seq(0, 20, by = 0.01)
js_c <- seq(-10, 10, by = 0.01)
fl_c <- fl(jl_c, tl, a, w)
fs_14_c <- fs_14(js_c, ts, a, w)
fs_17_c <- fs_17(jl_c, ts, a, w)
df_c <- data.frame(jl = jl_c,
                   js = js_c,
                   fl = fl_c,
                   fs_14 = fs_14_c,
                   fs_17 = fs_17_c,
                   eps = 0.12,
                   color1 = 1,
                   color3 = 3)

jl_d <- seq(0, 20, by = 1)
js_d <- seq(-10, 10, by = 1)
fl_d <- fl(jl_d, tl, a, w)
fs_14_d <- fs_14(js_d, ts, a, w)
fs_17_d <- fs_17(jl_d, ts, a, w)

df_d <- data.frame(jl = jl_d,
                   js = js_d,
                   fl = fl_d,
                   fs_14 = fs_14_d,
                   fs_17 = fs_17_d,
                   color2 = 2,
                   color4 = 4)

df_eps <- data.frame(x = 13.5,
                  y = min(fl_c):max(fl_c),
                  color3 = 3)


# Figure 1
########## value of terms plot, large-time
fig_1 <- ggplot() +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_line(data = df_c, size = 1,
            aes(x = jl, y = fl, color = factor(color1, levels = 1:4))) +
  geom_point(data = df_d[df_d[["jl"]] < 13.5,], size = 4, shape = 16,
             aes(x = jl, y = fl, color = factor(color2, levels = 1:4))) +
  geom_line(data = df_eps, linetype = "dashed", size = 1,
            aes(x = x, y = y, color = factor(color3, levels = 1:4))) +
  geom_point(data = df_d[df_d[["jl"]] > 13.5,], size = 4, shape = 17,
             aes(x = jl, y = fl, color = factor(color4, levels = 1:4))) +
  labs(x = "Index of the term",
       y = "Value of the term") +
  scale_color_manual(values = c("#decbc9", "#cb5041", "blue", "blue"),
                     name = NULL,
                     breaks = 1:4,
                     labels = c("continuous solution", "discretized terms",
                                "truncation index", "unnecessary terms")) +
  guides(color = guide_legend(override.aes = list(
    shape = c(NA, 16, NA, 17),
    linetype = c("solid", "blank", "dashed", "blank")))) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.border = element_blank(),
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
ggsave(paste0(img_dir, "terms_large.png"),
       plot = fig_1, width = 16, height = 9)



# Figure 2
########## value of terms plot, small-time, S_14 style
fig_2 <- ggplot() +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_line(data = df_c, size = 1,
            aes(x = js, y = fs_14, color = factor(color1, levels = 1:4))) +
  geom_point(data = df_d[abs(df_d[["fs_14"]]) > eps,], size = 4, shape = 16,
             aes(x = js, y = fs_14, color = factor(color2, levels = 1:4))) +
  geom_line(data = df_c, linetype = "dashed", size = 1,
            aes(x = js, y = eps, color = factor(color3, levels = 1:4))) +
  geom_line(data = df_c, linetype = "dashed", size = 1,
            aes(x = js, y = -eps, color = factor(color3, levels = 1:4))) +
  geom_point(data = df_d[abs(df_d[["fs_14"]]) < eps,], size = 4, shape = 17,
             aes(x = js, y = fs_14, color = factor(color4, levels = 1:4))) +
  labs(x = "Index of the term",
       y = "Value of the term") +
  scale_color_manual(values = c("#decbc9", "#cb5041", "blue", "blue"),
                     name = NULL,
                     breaks = 1:4,
                     labels = c("continuous solution", "discretized terms",
                                "error tolerance", "unnecessary terms")) +
  guides(color = guide_legend(override.aes = list(
    shape = c(NA, 16, NA, 17),
    linetype = c("solid", "blank", "dashed", "blank")))) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.border = element_blank(),
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
ggsave(paste0(img_dir, "terms_small_14.png"),
       plot = fig_2, width = 16, height = 9)


# Figure 3
########## value of terms plot, small-time, S_17 style
fig_3 <- ggplot() +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_point(data = df_d[abs(df_d[["fs_17"]]) > eps,], size = 4, shape = 16,
             aes(x = jl, y = fs_17, color = factor(color2, levels = 2:4))) +
  geom_line(data = df_c, linetype = "dashed", size = 1,
            aes(x = jl, y = eps, color = factor(color3, levels = 2:4))) +
  geom_line(data = df_c, linetype = "dashed", size = 1,
            aes(x = jl, y = -eps, color = factor(color3, levels = 2:4))) +
  geom_point(data = df_d[abs(df_d[["fs_17"]]) < eps,], size = 4, shape = 17,
             aes(x = jl, y = fs_17, color = factor(color4, levels = 2:4))) +
  labs(x = "Index of the term",
       y = "Value of the term") +
  scale_color_manual(values = c("#cb5041", "blue", "blue"),
                     name = NULL,
                     breaks = 2:4,
                     labels = c("discretized terms", "error tolerance",
                                "unnecessary terms")) +
  guides(color = guide_legend(override.aes = list(
    shape = c(16, NA, 17),
    linetype = c("blank", "dashed", "blank")))) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.border = element_blank(),
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
ggsave(paste0(img_dir, "terms_small_17.png"),
       plot = fig_3, width = 16, height = 9)
