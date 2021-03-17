# Figure S1

# This file produces the supplementary results (and plots) that pertain to the
# implementations that combine the SWSE small-time method with the
# Navarro et al. large-time method on a predefined parameter space, as shown in
# Section 2.1 of the supplementary analysis document of the fddm paper:
# "Determining the Default Behavior of our Heuristic Switching Mechanism, delta"

library("fddm")
library("microbenchmark")
library("reshape2")
library("ggplot2")
save_dir <- "paper_analysis/supplemental_analysis/2_delta/results/"
img_dir <- "paper_analysis/images/supplemental_analysis/2_delta/"



##### Benchmark Function #######################################################
rt_benchmark_vec <- function(RT, resp, V, A, t0 = 1e-4, W = 0.5, SV = 0.0,
                             err_tol = 1e-6, times = 100, unit = "ns") {

  fnames <- c("fb_SWSE_17_0", "fb_SWSE_14_0", "fb_SWSE_17_1", "fb_SWSE_14_1",
              "fb_SWSE_17_2", "fb_SWSE_14_2", "fb_SWSE_17_3", "fb_SWSE_14_3",
              "fb_SWSE_17_4", "fb_SWSE_14_4", "fb_SWSE_17_5", "fb_SWSE_14_5",
              "fb_SWSE_17_6", "fb_SWSE_14_6", "fb_SWSE_17_7", "fb_SWSE_14_7")
  nf <- length(fnames) # number of functions being benchmarked
  nV <- length(V)
  nA <- length(A)
  nW <- length(W)
  nSV <- length(SV)
  resp <- rep(resp, length(RT)) # for RWiener

  # Initialize the dataframe to contain the microbenchmark results
  mbm_res <- data.frame(matrix(ncol = 4+nf, nrow = nV*nA*nW*nSV))
  colnames(mbm_res) <- c('V', 'A', 'W', 'SV', fnames)
  row_idx <- 1

  # Loop through each combination of parameters and record microbenchmark results
  for (v in 1:nV) {
    for (a in 1:nA) {
      for (w in 1:nW) {
        for (sv in 1:nSV) {
          mbm <- microbenchmark(
            fb_SWSE_17_0 = dfddm(rt = RT, response = resp, a = A[a],
                                 v = V[v], t0 = t0, w = W[w],
                                 log = FALSE, n_terms_small = "SWSE",
                                 summation_small = "2017", scale = "both",
                                 max_terms_large = 0, err_tol = err_tol),
            fb_SWSE_14_0 = dfddm(rt = RT, response = resp, a = A[a],
                                 v = V[v], t0 = t0, w = W[w],
                                 log = FALSE, n_terms_small = "SWSE",
                                 summation_small = "2014", scale = "both",
                                 max_terms_large = 0, err_tol = err_tol),
            fb_SWSE_17_1 = dfddm(rt = RT, response = resp, a = A[a],
                                 v = V[v], t0 = t0, w = W[w],
                                 log = FALSE, n_terms_small = "SWSE",
                                 summation_small = "2017", scale = "both",
                                 max_terms_large = 1, err_tol = err_tol),
            fb_SWSE_14_1 = dfddm(rt = RT, response = resp, a = A[a],
                                 v = V[v], t0 = t0, w = W[w],
                                 log = FALSE, n_terms_small = "SWSE",
                                 summation_small = "2014", scale = "both",
                                 max_terms_large = 1, err_tol = err_tol),
            fb_SWSE_17_2 = dfddm(rt = RT, response = resp, a = A[a],
                                 v = V[v], t0 = t0, w = W[w],
                                 log = FALSE, n_terms_small = "SWSE",
                                 summation_small = "2017", scale = "both",
                                 max_terms_large = 2, err_tol = err_tol),
            fb_SWSE_14_2 = dfddm(rt = RT, response = resp, a = A[a],
                                 v = V[v], t0 = t0, w = W[w],
                                 log = FALSE, n_terms_small = "SWSE",
                                 summation_small = "2014", scale = "both",
                                 max_terms_large = 2, err_tol = err_tol),
            fb_SWSE_17_3 = dfddm(rt = RT, response = resp, a = A[a],
                                 v = V[v], t0 = t0, w = W[w],
                                 log = FALSE, n_terms_small = "SWSE",
                                 summation_small = "2017", scale = "both",
                                 max_terms_large = 3, err_tol = err_tol),
            fb_SWSE_14_3 = dfddm(rt = RT, response = resp, a = A[a],
                                 v = V[v], t0 = t0, w = W[w],
                                 log = FALSE, n_terms_small = "SWSE",
                                 summation_small = "2014", scale = "both",
                                 max_terms_large = 3, err_tol = err_tol),
            fb_SWSE_17_4 = dfddm(rt = RT, response = resp, a = A[a],
                                 v = V[v], t0 = t0, w = W[w],
                                 log = FALSE, n_terms_small = "SWSE",
                                 summation_small = "2017", scale = "both",
                                 max_terms_large = 4, err_tol = err_tol),
            fb_SWSE_14_4 = dfddm(rt = RT, response = resp, a = A[a],
                                 v = V[v], t0 = t0, w = W[w],
                                 log = FALSE, n_terms_small = "SWSE",
                                 summation_small = "2014", scale = "both",
                                 max_terms_large = 4, err_tol = err_tol),
            fb_SWSE_17_5 = dfddm(rt = RT, response = resp, a = A[a],
                                 v = V[v], t0 = t0, w = W[w],
                                 log = FALSE, n_terms_small = "SWSE",
                                 summation_small = "2017", scale = "both",
                                 max_terms_large = 5, err_tol = err_tol),
            fb_SWSE_14_5 = dfddm(rt = RT, response = resp, a = A[a],
                                 v = V[v], t0 = t0, w = W[w],
                                 log = FALSE, n_terms_small = "SWSE",
                                 summation_small = "2014", scale = "both",
                                 max_terms_large = 5, err_tol = err_tol),
            fb_SWSE_17_6 = dfddm(rt = RT, response = resp, a = A[a],
                                 v = V[v], t0 = t0, w = W[w],
                                 log = FALSE, n_terms_small = "SWSE",
                                 summation_small = "2017", scale = "both",
                                 max_terms_large = 6, err_tol = err_tol),
            fb_SWSE_14_6 = dfddm(rt = RT, response = resp, a = A[a],
                                 v = V[v], t0 = t0, w = W[w],
                                 log = FALSE, n_terms_small = "SWSE",
                                 summation_small = "2014", scale = "both",
                                 max_terms_large = 6, err_tol = err_tol),
            fb_SWSE_17_7 = dfddm(rt = RT, response = resp, a = A[a],
                                 v = V[v], t0 = t0, w = W[w],
                                 log = FALSE, n_terms_small = "SWSE",
                                 summation_small = "2017", scale = "both",
                                 max_terms_large = 7, err_tol = err_tol),
            fb_SWSE_14_7 = dfddm(rt = RT, response = resp, a = A[a],
                                 v = V[v], t0 = t0, w = W[w],
                                 log = FALSE, n_terms_small = "SWSE",
                                 summation_small = "2014", scale = "both",
                                 max_terms_large = 7, err_tol = err_tol),
            times = times, unit = unit)
          # add the v, a, and w values to the dataframe
          mbm_res[row_idx, 1] <- V[v]
          mbm_res[row_idx, 2] <- A[a]
          mbm_res[row_idx, 3] <- W[w]
          mbm_res[row_idx, 4] <- SV[sv]
          # add the median microbenchmark results to the dataframe
          for (i in 1:nf) {
            mbm_res[row_idx, 4+i] <- median(mbm[mbm[,1] == fnames[i],2])
          }
          # iterate start value
          row_idx = row_idx + 1
        }
      }
    }
  }
  return(mbm_res)
}




##### Parameter Space 2 (RT:0-2) Benchmorks ####################################
RT <- seq(0.1, 2, by = 0.1)
A <- seq(0.5, 3.5, by = 0.5)
V <- c(-5, -2, 0, 2, 5)
t0 <- 1e-4 # must be nonzero for RWiener
W <- seq(0.3, 0.7, by = 0.1)
SV <- c(0, 1, 2, 3.5)
err_tol <- 1e-6 # this is the setting from rtdists

bm_vec_2 <- rt_benchmark_vec(RT = RT, resp = "lower", V = V, A = A, t0 = t0,
                             W = W, SV = SV, err_tol = err_tol,
                             times = 1000, unit = "ns")

save(bm_vec_2, compress = "xz", compression_level = 9,
     file = paste0(save_dir, "delta_bm_0-2.Rds"))

# uncomment the following line if loading pre-run benchmark data,
# will load into variable named "bm_vec_2"
# load(paste0(save_dir, "delta_bm_0-2.Rds"))

# Plot Results
t_idx <- match("SV", colnames(bm_vec_2))
bm_vec_2[, -seq_len(t_idx)] <- bm_vec_2[, -seq_len(t_idx)]/1000 # convert to microseconds
mbm_vec_2 <- melt(bm_vec_2, measure.vars = -seq_len(t_idx),
                variable.name = "FuncName", value.name = "time")

Names_vec <- c("fb_SWSE_17_0", "fb_SWSE_14_0", "fb_SWSE_17_1", "fb_SWSE_14_1",
               "fb_SWSE_17_2", "fb_SWSE_14_2", "fb_SWSE_17_3", "fb_SWSE_14_3",
               "fb_SWSE_17_4", "fb_SWSE_14_4", "fb_SWSE_17_5", "fb_SWSE_14_5",
               "fb_SWSE_17_6", "fb_SWSE_14_6", "fb_SWSE_17_7", "fb_SWSE_14_7")
Color_vec <- c("#e000b4", "#ff99eb", "#e68a00", "#ffb366",
               "#006699", "#66ccff", "#9900cc", "#cc99ff",
               "#c2a500", "#d7db42", "#336600", "#33cc33",
               "#996633", "#ff9999", "#ff5050", "#990000")

# Figure S1
fig_s1 <- ggplot(data = mbm_vec_2,
                 aes(x = factor(FuncName, levels = Names_vec), y = time,
                     color = factor(FuncName, levels = Names_vec),
                     fill = factor(FuncName, levels = Names_vec))) +
  geom_violin(trim = TRUE, alpha = 0.5) +
  scale_x_discrete(labels = c(bquote("0, " ~ S[17]), bquote("0, " ~ S[14]),
                              bquote("1, " ~ S[17]), bquote("1, " ~ S[14]),
                              bquote("2, " ~ S[17]), bquote("2, " ~ S[14]),
                              bquote("3, " ~ S[17]), bquote("3, " ~ S[14]),
                              bquote("4, " ~ S[17]), bquote("4, " ~ S[14]),
                              bquote("5, " ~ S[17]), bquote("5, " ~ S[14]),
                              bquote("6, " ~ S[17]), bquote("6, " ~ S[14]),
                              bquote("7, " ~ S[17]), bquote("7, " ~ S[14]))) +
  scale_color_manual(values = Color_vec, guide = FALSE) +
  scale_fill_manual(values = Color_vec, guide = FALSE) +
  geom_boxplot(width = 0.15, fill = "white", alpha = 0.5) +
  stat_summary(fun = mean, geom = "errorbar",
               aes(ymax = ..y.., ymin = ..y..),
               width = .35, linetype = "dashed") +
  labs(x = bquote(delta ~ ", and Summation Style"), y = "Time (microseconds)") +
  theme_bw() +
  theme(panel.border = element_blank(),
        axis.text.x = element_text(size = 16, angle = 90,
                                   vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 20,
                                    margin = margin(10, 5, 5, 5)),
        axis.title.y = element_text(size = 20),
        legend.position = "none")
ggsave(paste0(img_dir, "delta_bm_vec_0-2.png"),
       plot = fig_s1, width = 16, height = 9)
