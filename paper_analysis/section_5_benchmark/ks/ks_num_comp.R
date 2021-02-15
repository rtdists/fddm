library("Rcpp")
library("reshape2")
library("ggplot2")
save_dir <- "paper_analysis/section_5_benchmark/ks/"

# Load the ks functions (from C++, using Rcpp)
# the functions are:
  # ks_Gon(t, w, eps)
  # ks_Nav(t, w, eps) (NOTE: ks_Nav does not use the input `w`)
	# ks_SWSE_14(t, w, eps)
	# ks_SWSE_17(t, w, eps)
sourceCpp(paste0(save_dir, "ks_vec.cpp"))

ks_SWSE_14(0.9, 0, 1e-6)
ks_SWSE_17(0.9, 0, 1e-6)
ks_Gon(0.9, 0, 1e-6)


##### Testing Function #########################################################
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
		ks_SWSE_14 = ks_SWSE_14(t, w, eps),
		ks_SWSE_17 = ks_SWSE_17(t, w, eps)
  )

  return(res)
}



##### Run and Save the Results #################################################
t <- seq(0.1, 10, by = 0.1)
w <- seq(0, 1, by = 0.1)
eps <- 10^(-1:-10)

ks_num <- ks_num_test(t, w, eps)
save(ks_num, compress = "xz", compression_level = 9,
     file = paste0(save_dir, "ks_num.Rds"))

tt <- 1.3
ww <- 0.1
eeps <- 0.1
gamma <- -1 / (2 * tt);
ww * exp(gamma * ww*ww)
sqrt(tt)/2 - ww/2
sqrt(tt) - ww

ks_num[ks_num[["ks_SWSE_17"]] - ks_num[["ks_SWSE_14"]] > 0, ]
hist(ks_num[ks_num[["ks_SWSE_17"]] - ks_num[["ks_SWSE_14"]] > 0, "w"])
all(ks_num[ks_num[["ks_SWSE_17"]] - ks_num[["ks_SWSE_14"]] > 0, "ks_SWSE_14"] == 1)

##### Plot Heatmap of Results ##################################################
# uncomment the following line if loading pre-run benchmark data,
# will load into variable named "ks_num"
# load(paste0(save_dir, "ks_num.Rds"))

ks_num[["ks_CGN"]] <- pmin(ks_num[["ks_Gondan_14"]], ks_num[["ks_Navarro_09"]])
ks_num[["ks_diff_GN"]] <- ks_num[["ks_Gondan_14"]] - ks_num[["ks_Navarro_09"]]
ks_num[["ks_diff_SC"]] <- ks_num[["ks_SWSE_17"]] - ks_num[["ks_CGN"]]


sum(ks_num[["ks_diff_GN"]])
sum(ks_num[["ks_diff_SC"]])

w <- sort(unique(ks_num[["w"]]))

# Gondan - Navarro
ks_diffs_GN <- sort(unique(ks_num[["ks_diff_GN"]]))
lbl_GN <- as.character(ks_diffs_GN)
lbl_GN[ks_diffs_GN > 0] <- paste0("+", lbl_GN[ks_diffs_GN > 0])

for (i in 1:length(w)) {
  ggplot(ks_num[ks_num[["w"]] == w[i], ]) +
    geom_tile(aes(x = t, y = eps, fill = ks_diff_GN)) +
    scale_fill_gradient2(low = "#c2a500",
                         mid = "#ffffff",
                         high = "#336600",
                         midpoint = 0,
                         name = bquote(k[s] ~ "diff"),
                         breaks = ks_diffs_GN,
                         labels = lbl_GN) +
    scale_y_log10() +
    labs(title = bquote("Differences in " ~ k[s] ~ " calculation for Gondan (2014) - Navarro (2009)"),
         subtitle = paste0("w = ", w[i], "\n",
                           "Cumulative difference = ", sum(ks_num[ks_num[["w"]] == w[i], "ks_diff_GN"])),
         x = bquote(frac(rt, a^2) ~ ", effective response time"),
         y = bquote(epsilon ~ ", desired precision, " ~ log[10])) +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          plot.title = element_text(size = 20),
          plot.subtitle = element_text(size = 17),
          axis.text.x = element_text(size = 16),
          axis.text.y = element_text(size = 16),
          axis.title.x = element_text(size = 18),
          axis.title.y = element_text(size = 18),
          strip.text = element_text(size = 14),
          strip.background = element_rect(fill = "white"))
  ggsave(paste0("paper_analysis/images/ks_heatmap/Gondan_Navarro/w_", i, ".png"),
         width = 16, height = 9)
}


# SWSE - min(Gondan, Navarro)
ks_diffs_SC <- sort(unique(ks_num[["ks_diff_SC"]]))
lbl_SC <- as.character(ks_diffs_SC)
lbl_SC[ks_diffs_SC > 0] <- paste0("+", lbl_SC[ks_diffs_SC > 0])

for (i in 1:length(w)) {
  ggplot(ks_num[ks_num[["w"]] == w[i], ]) +
    geom_tile(aes(x = t, y = eps, fill = ks_diff_SC)) +
    scale_fill_gradient2(low = "#9900cc",
                         mid = "#ffffff",
                         high = "#336600",
                         midpoint = 0,
                         name = bquote(k[s] ~ "diff"),
                         breaks = ks_diffs_SC,
                         labels = lbl_SC) +
    scale_y_log10() +
    labs(title = bquote("Differences in " ~ k[s] ~ " calculation for SWSE - min(Gondan (2014), Navarro (2009) )"),
         subtitle = paste0("w = ", w[i], "\n",
                           "Cumulative difference = ", sum(ks_num[ks_num[["w"]] == w[i], "ks_diff_SC"])),
         x = bquote(frac(rt, a^2) ~ ", effective response time"),
         y = bquote(epsilon ~ ", desired precision, " ~ log[10])) +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          plot.title = element_text(size = 20),
          plot.subtitle = element_text(size = 17),
          axis.text.x = element_text(size = 16),
          axis.text.y = element_text(size = 16),
          axis.title.x = element_text(size = 18),
          axis.title.y = element_text(size = 18),
          strip.text = element_text(size = 14),
          strip.background = element_rect(fill = "white"))
  ggsave(paste0("paper_analysis/images/ks_heatmap/SWSE_CGN/w_", i, ".png"),
         width = 16, height = 9)
}
