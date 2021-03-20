# Figure 7

# This file produces the results (and plot) that pertain to the differences in
# speed between our implementations of the small-time approximation methods
# provided by Navarro et al. and Gondan et al., as shown in
# Section 5.3.1 of the fddm paper:
# "Benchmarking All Implementations/Predefined Parameter Spaces"

library("Rcpp")
library("reshape2")
library("ggplot2")
save_dir <- "paper_analysis/sec_5_benchmark/3_all_implementations/"
img_dir <- "paper_analysis/images/sec_5_benchmark/"

# Load the ks functions (from C++, using Rcpp)
# the functions are:
  # ks_Gon(t, w, eps)
  # ks_Nav(t, w, eps) (NOTE: ks_Nav does not use the input `w`)
sourceCpp(paste0(save_dir, "ks.cpp"))



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
    ks_Navarro_09 = ks_Nav(t, w, eps),
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

# Figure 7
fig_7 <- ggplot(ks_num[ks_num[["w"]] == 0.5, ]) +
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
  labs(subtitle = paste0("Cumulative difference = ",
                         sum(ks_num[ks_num[["w"]] == 0.5, "ks_diff_GN"])),
       x = bquote(frac(rt, a^2) ~ ", effective response time"),
       y = bquote(epsilon ~ ", desired precision, " ~ log[10])) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        plot.subtitle = element_text(size = 17),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        strip.text = element_text(size = 14),
        strip.background = element_rect(fill = "white"))
ggsave(paste0(img_dir, "ks_Gon_Nav_w_5.png"),
       plot = fig_7, width = 16, height = 9)
