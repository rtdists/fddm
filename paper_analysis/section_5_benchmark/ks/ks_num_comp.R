library("Rcpp")
library("reshape2")
library("ggplot2")
save_dir <- "paper_analysis/section_5_benchmark/ks/"

# Load the ks functions (from C++, using Rcpp)
# the functions are:
  # ks_Gon(t, w, eps)
  # ks_Nav(t, w, eps) (NOTE: ks_Nav does not use the input `w`)
sourceCpp(paste0(save_dir, "ks_vec.cpp"))





##### Testing Function #########################################################
t <- seq(0.1, 10, by = 0.1)
w <- seq(0, 1, by = 0.1)
eps <- 10^(-1:-10)




ks_num_test <- function(t, w, eps) {
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



for (ti in 0:(nt-1)) {
  for (wi in 0:(nw-1)) {
    for (ei in 0:(ne-1)) {
      i <- ti*(nw*ne) + wi*ne + ei
      cat(i, " = ", ti*(nw*ne), " + ", wi*ne, " + ", ei, "\n")
    }
  }
}



##### Run and Save the Results #################################################
ks_num <- ks_num_test(t, w, eps)
save(ks_num, compress = "xz", compression_level = 9,
     file = paste0(save_dir, "ks_num.Rds"))



##### Plot Heatmap of Results ##################################################
# uncomment the following line if loading pre-run benchmark data,
# will load into variable named "ks_num"
# load(paste0(save_dir, "ks_num.Rds"))

ks_num[["ks_diff"]] <- ks_num[["ks_Gondan_14"]] - ks_num[["ks_Navarro_09"]]
ks_diffs <- sort(unique(ks_num[["ks_diff"]]))
lbl <- as.character(ks_diffs)
lbl[ks_diffs > 0] <- paste0("+", lbl[ks_diffs > 0])
w <- sort(unique(ks_num[["w"]]))

for (i in 1:length(w)) {
  ggplot(ks_num[ks_num[["w"]] == w[i], ]) +
    geom_tile(aes(x = t, y = eps, fill = ks_diff)) +
    scale_fill_gradient2(low = "#990000",
                         mid = "#ffffff",
                         high = "#000099",
                         midpoint = 0,
                         name = bquote(k[s] ~ "diff"),
                         breaks = ks_diffs,
                         labels = lbl) +
    scale_x_continuous(breaks = unique(ks_num[["t"]]), labels = xlab) +
    scale_y_log10() +
    labs(title = bquote("Median microbenchmark results for " ~ k[s] ~ " calculation"),
         x = bquote(frac(rt, a^2) ~ ", effective response time"),
         y = bquote(epsilon ~ ", desired precision")) +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          plot.title = element_text(size = 20,
                                    margin = margin(5, 5, 30, 5, "pt")),
          axis.text.x = element_text(size = 16),
          axis.text.y = element_text(size = 16),
          axis.title.x = element_text(size = 18),
          axis.title.y = element_text(size = 18),
          strip.text = element_text(size = 14),
          strip.background = element_rect(fill = "white"))
  ggsave(paste0("paper_analysis/images/ks_heatmap/w_", i, ".png"),
         width = 16, height = 9)
}
