# Figures S6, S7

# This file produces the supplementary results (and plot) that pertain to the
# differences in speed between our implementations of the small-time
# approximation methods given by Navarro et al. and Gondan et al., as shown in
# Section 2.2 of the supplementary analysis document of the fddm paper:
# "Benchmarking All Implementations/Predefined Parameter Spaces"

library("Rcpp")
library("microbenchmark")
library("reshape2")
library("ggplot2")
save_dir <- "paper_analysis/supplemental_analysis/3_all_implementations/"
img_dir <- "paper_analysis/images/supplemental_analysis/3_all_implementations/"

# Load the ks functions (from C++, using Rcpp)
# the functions are:
  # ks_Gon(t, w, eps)
  # ks_Nav(t, w, eps) (NOTE: ks_Nav does not use the input `w`)
sourceCpp(paste0(save_dir, "ks.cpp"))



########## Comparing ks Results ################################################
### ks Sampling Function
ks_num_test <- function(t, w, eps) {
  t <- t[t > 0]
  nt <- length(t)
  nw <- length(w)
  ne <- length(eps)

  # Initialize the dataframe and fill with the results
  res <- data.frame(
    t = rep(t, each = nw * ne),
    w = rep(w, each = ne, times = nt),
    eps = rep(eps, each = 1, times = nw * nt),
    ks_Gondan_14 = ks_Gon(t, w, eps),
    ks_Navarro_09 = ks_Nav(t, w, eps)
  )

  return(res)
}


### Run and save the results (doesn't take long)
t <- seq(0.1, 10, by = 0.1)
w <- seq(0.1, 0.9, by = 0.1)
eps <- 10^(-2:-10)

ks_num <- ks_num_test(t, w, eps)
save(ks_num, compress = "xz", compression_level = 9,
     file = paste0(save_dir, "results/ks_num.Rds"))


### Plot heatmap of results
# uncomment the following line if loading pre-run benchmark data,
# will load into variable named "ks_num"
# load(paste0(save_dir, "results/ks_num.Rds"))

ks_num[["ks_diff_GN"]] <- ks_num[["ks_Gondan_14"]] - ks_num[["ks_Navarro_09"]]

max_diff_GN <- max(abs(ks_num[["ks_diff_GN"]]))
breaks_GN <- trunc(seq(from = -max_diff_GN, to = max_diff_GN, length.out = 5))
lbl_GN <- as.character(breaks_GN)
lbl_GN[breaks_GN > 0] <- paste0("+", lbl_GN[breaks_GN > 0])

w <- c(0.1, 0.3, 0.7, 0.9)

# Figure 7
for (i in 1:length(w)) {
  fig_7 <- ggplot(ks_num[ks_num[["w"]] >= w[i]-0.01 & ks_num[["w"]] <= w[i]+0.01, ]) +
    geom_tile(aes(x = t, y = eps, fill = ks_diff_GN)) +
    scale_fill_gradient2(low = "#b34d4d",
                         mid = "#ffffff",
                         high = "#4d80b3",
                         midpoint = 0,
                         name = bquote(k[s] ~ "difference"),
                         limits = c(-max_diff_GN, max_diff_GN),
                         breaks = breaks_GN,
                         labels = lbl_GN) +
    scale_y_log10() +
    labs(subtitle = paste0("w = ", w[i], ";\n",
                           "Cumulative difference = ",
                           sum(ks_num[ks_num[["w"]] >= w[i]-0.01 & ks_num[["w"]] <= w[i]+0.01, "ks_diff_GN"])),
         x = bquote(frac(t, a^2) ~ ", effective response time"),
         y = bquote(epsilon ~ ", desired precision, " ~ log[10])) +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          plot.subtitle = element_text(size = 18),
          axis.text.x = element_text(size = 16),
          axis.text.y = element_text(size = 16),
          axis.title.x = element_text(size = 20,
                                      margin = margin(10, 0, 0, 0)),
          axis.title.y = element_text(size = 20,
                                      margin = margin(0, 10, 0, 0)),
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 13))
  ggsave(paste0(img_dir, "ks_heatmap/ks_Gon_Nav_w_", i, ".png"),
         plot = fig_7, width = 16, height = 9)
}



########## Comparing ks Benchmarks #############################################
### Benchmark Function
ks_benchmark <- function(t_range, dt, def = 1000, w, eps,
                         times = 1000, unit = "ns") {
  fnames <- c("ks_Gondan_14", "ks_Navarro_09")
  nf <- length(fnames) # number of functions being benchmarked

  if (length(w) != 1) {stop("length of w needs to be 1")}
  if (length(eps) != 1) {stop("length of eps needs to be 1")}

  t_min <- min(t_range)
  if (t_min <= 0) {t_min <- dt}
  t_max <- max(t_range)
  nt <- ceiling((t_max - t_min) / dt)

  # Initialize the dataframe to contain the microbenchmark results
  mbm_res <- data.frame(matrix(ncol = 3+nf, nrow = nt))
  colnames(mbm_res) <- c('t_start', 't_mid', 't_stop', fnames)

  # Loop through each combination of parameters
  for (i in 1:nt) {
    # slice t
    t_start <- t_min + (i-1)*dt
    t_stop <- t_start + dt
    t <- seq(from = t_start, to = t_stop, length.out = def)
    mbm <- microbenchmark(
      # record microbenchmark results
      ks_Gondan_14 = ks_Gon(t, w, eps),
      ks_Navarro_09 = ks_Nav(t, w, eps),
      times = times, unit = unit)
    # add the t, w, and eps values to the dataframe
    mbm_res[i, 1] <- t_start
    mbm_res[i, 2] <- t_start + (dt / 2)
    mbm_res[i, 3] <- t_stop
    # add the median microbenchmark results to the dataframe
    for (fi in 1:nf) {
      mbm_res[i, 3+fi] <- median(mbm[mbm[,1] == fnames[fi], 2])
    }
  }
  return(mbm_res)
}


### Run and Save the Benchmark Results
ks_bm <- ks_benchmark(t_range = c(0, 10), dt = 0.02, def = 1000,
                      w = 0.5, eps = 1e-6, times = 1000)
save(ks_bm, compress = "xz", compression_level = 9,
     file = paste0(save_dir, "ks_bm.Rds"))


### Plot Benchmark Timing Results
# Figure S7

Names <- c("ks_Gondan_14", "ks_Navarro_09")
Color <- c("#b34d4d", "#4d80b3")

# uncomment the following line if loading pre-run benchmark data,
# will load into variable named "ks_bm"
# load(paste0(save_dir, "results/ks_bm-windows.Rds"))

t_idx <- match("t_stop", colnames(ks_bm))
mbm_ks <- melt(ks_bm, measure.vars = -seq_len(t_idx),
               variable.name = "FuncName", value.name = "time")
mbm_ks[["time"]] <- mbm_ks[["time"]] * 1e-3 # convert nano to microseconds

fig_s7a <- ggplot(mbm_ks, aes(x = t_mid, y = time,
                              color = factor(FuncName, levels = Names))) +
  geom_line(size = 1) +
  scale_color_manual(values = Color,
                     labels = c(bquote(k[s] ~ "Gondan (2014)"),
                                bquote(k[s] ~ "Navarro (2009)")),
                     name = "Method:") +
  guides(color = guide_legend(override.aes = list(size = c(2,2)))) +
  labs(title = "", # need title because of position of legend
       x = bquote(frac(t, a^2) ~ ", effective response time"),
       y = "Time (microseconds)") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        plot.title = element_text(size = 20,
                                  margin = margin(0, 0, 0, 0)),
        axis.text.x = element_text(size = 16, angle = 90,
                                   vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.position = c(1, 1),
        legend.justification = c(1, 0),
        legend.box = "horizontal",
        legend.direction = "horizontal",
        legend.background = element_rect(fill = "transparent"),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 13))
ggsave(paste0(img_dir, "ks_med_t-windows.png"),
       plot = fig_s7a, width = 16, height = 9)

# uncomment the following line if loading pre-run benchmark data,
# will load into variable named "ks_bm"
# load(paste0(save_dir, "ks_bm-linux.Rds"))

t_idx <- match("t_stop", colnames(ks_bm))
mbm_ks <- melt(ks_bm, measure.vars = -seq_len(t_idx),
               variable.name = "FuncName", value.name = "time")
mbm_ks[["time"]] <- mbm_ks[["time"]] * 1e-3 # convert nano to microseconds

fig_s7b <- ggplot(mbm_ks, aes(x = t_mid, y = time,
                              color = factor(FuncName, levels = Names))) +
  geom_line(size = 1) +
  scale_color_manual(values = Color,
                     labels = c(bquote(k[s] ~ "Gondan (2014)"),
                                bquote(k[s] ~ "Navarro (2009)")),
                     name = "Method:") +
  guides(color = guide_legend(override.aes = list(size = c(2,2)))) +
  labs(title = "", # need title because of position of legend
       x = bquote(frac(t, a^2) ~ ", effective response time"),
       y = "Time (microseconds)") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        plot.title = element_text(size = 20,
                                  margin = margin(0, 0, 0, 0)),
        axis.text.x = element_text(size = 16, angle = 90,
                                   vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.position = c(1, 1),
        legend.justification = c(1, 0),
        legend.box = "horizontal",
        legend.direction = "horizontal",
        legend.background = element_rect(fill = "transparent"),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 13))
ggsave(paste0(img_dir, "ks_med_t-linux.png"),
       plot = fig_s7b, width = 16, height = 9)
