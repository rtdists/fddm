library("fddm")
library("microbenchmark")
library("reshape2")
library("ggplot2")
save_dir <- "paper_analysis/extra_analysis/"



##### Benchmark Functions ######################################################
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

rt_benchmark_ind <- function(RT, resp, V, A, t0 = 1e-4, W = 0.5, SV = 0.0,
                             err_tol = 1e-6, times = 100, unit = "ns") {
  fnames <- c("fb_SWSE_17_0", "fb_SWSE_14_0", "fb_SWSE_17_1", "fb_SWSE_14_1",
              "fb_SWSE_17_2", "fb_SWSE_14_2", "fb_SWSE_17_3", "fb_SWSE_14_3",
              "fb_SWSE_17_4", "fb_SWSE_14_4", "fb_SWSE_17_5", "fb_SWSE_14_5",
              "fb_SWSE_17_6", "fb_SWSE_14_6", "fb_SWSE_17_7", "fb_SWSE_14_7")
  nf <- length(fnames) # number of functions being benchmarked
  nRT <- length(RT)
  nV <- length(V)
  nA <- length(A)
  nW <- length(W)
  nSV <- length(SV)

  # Initialize the dataframe to contain the microbenchmark results
  mbm_res <- data.frame(matrix(ncol = 5+nf, nrow = nRT*nV*nA*nW*nSV))
  colnames(mbm_res) <- c('RT', 'V', 'A', 'W', 'SV', fnames)
  row_idx <- 1

  # Loop through each combination of parameters and record microbenchmark results
  for (rt in 1:nRT) {
    for (v in 1:nV) {
      for (a in 1:nA) {
        for (w in 1:nW) {
          for (sv in 1:nSV) {
            mbm <- microbenchmark(
              fb_SWSE_17_0 = dfddm(rt = RT[rt], response = resp, a = A[a],
                                   v = V[v], t0 = t0, w = W[w],
                                   log = FALSE, n_terms_small = "SWSE",
                                   summation_small = "2017", scale = "both",
                                   max_terms_large = 0, err_tol = err_tol),
              fb_SWSE_14_0 = dfddm(rt = RT[rt], response = resp, a = A[a],
                                   v = V[v], t0 = t0, w = W[w],
                                   log = FALSE, n_terms_small = "SWSE",
                                   summation_small = "2014", scale = "both",
                                   max_terms_large = 0, err_tol = err_tol),
              fb_SWSE_17_1 = dfddm(rt = RT[rt], response = resp, a = A[a],
                                   v = V[v], t0 = t0, w = W[w],
                                   log = FALSE, n_terms_small = "SWSE",
                                   summation_small = "2017", scale = "both",
                                   max_terms_large = 1, err_tol = err_tol),
              fb_SWSE_14_1 = dfddm(rt = RT[rt], response = resp, a = A[a],
                                   v = V[v], t0 = t0, w = W[w],
                                   log = FALSE, n_terms_small = "SWSE",
                                   summation_small = "2014", scale = "both",
                                   max_terms_large = 1, err_tol = err_tol),
              fb_SWSE_17_2 = dfddm(rt = RT[rt], response = resp, a = A[a],
                                   v = V[v], t0 = t0, w = W[w],
                                   log = FALSE, n_terms_small = "SWSE",
                                   summation_small = "2017", scale = "both",
                                   max_terms_large = 2, err_tol = err_tol),
              fb_SWSE_14_2 = dfddm(rt = RT[rt], response = resp, a = A[a],
                                   v = V[v], t0 = t0, w = W[w],
                                   log = FALSE, n_terms_small = "SWSE",
                                   summation_small = "2014", scale = "both",
                                   max_terms_large = 2, err_tol = err_tol),
              fb_SWSE_17_3 = dfddm(rt = RT[rt], response = resp, a = A[a],
                                   v = V[v], t0 = t0, w = W[w],
                                   log = FALSE, n_terms_small = "SWSE",
                                   summation_small = "2017", scale = "both",
                                   max_terms_large = 3, err_tol = err_tol),
              fb_SWSE_14_3 = dfddm(rt = RT[rt], response = resp, a = A[a],
                                   v = V[v], t0 = t0, w = W[w],
                                   log = FALSE, n_terms_small = "SWSE",
                                   summation_small = "2014", scale = "both",
                                   max_terms_large = 3, err_tol = err_tol),
              fb_SWSE_17_4 = dfddm(rt = RT[rt], response = resp, a = A[a],
                                   v = V[v], t0 = t0, w = W[w],
                                   log = FALSE, n_terms_small = "SWSE",
                                   summation_small = "2017", scale = "both",
                                   max_terms_large = 4, err_tol = err_tol),
              fb_SWSE_14_4 = dfddm(rt = RT[rt], response = resp, a = A[a],
                                   v = V[v], t0 = t0, w = W[w],
                                   log = FALSE, n_terms_small = "SWSE",
                                   summation_small = "2014", scale = "both",
                                   max_terms_large = 4, err_tol = err_tol),
              fb_SWSE_17_5 = dfddm(rt = RT[rt], response = resp, a = A[a],
                                   v = V[v], t0 = t0, w = W[w],
                                   log = FALSE, n_terms_small = "SWSE",
                                   summation_small = "2017", scale = "both",
                                   max_terms_large = 5, err_tol = err_tol),
              fb_SWSE_14_5 = dfddm(rt = RT[rt], response = resp, a = A[a],
                                   v = V[v], t0 = t0, w = W[w],
                                   log = FALSE, n_terms_small = "SWSE",
                                   summation_small = "2014", scale = "both",
                                   max_terms_large = 5, err_tol = err_tol),
              fb_SWSE_17_6 = dfddm(rt = RT[rt], response = resp, a = A[a],
                                   v = V[v], t0 = t0, w = W[w],
                                   log = FALSE, n_terms_small = "SWSE",
                                   summation_small = "2017", scale = "both",
                                   max_terms_large = 6, err_tol = err_tol),
              fb_SWSE_14_6 = dfddm(rt = RT[rt], response = resp, a = A[a],
                                   v = V[v], t0 = t0, w = W[w],
                                   log = FALSE, n_terms_small = "SWSE",
                                   summation_small = "2014", scale = "both",
                                   max_terms_large = 6, err_tol = err_tol),
              fb_SWSE_17_7 = dfddm(rt = RT[rt], response = resp, a = A[a],
                                   v = V[v], t0 = t0, w = W[w],
                                   log = FALSE, n_terms_small = "SWSE",
                                   summation_small = "2017", scale = "both",
                                   max_terms_large = 7, err_tol = err_tol),
              fb_SWSE_14_7 = dfddm(rt = RT[rt], response = resp, a = A[a],
                                   v = V[v], t0 = t0, w = W[w],
                                   log = FALSE, n_terms_small = "SWSE",
                                   summation_small = "2014", scale = "both",
                                   max_terms_large = 7, err_tol = err_tol),
              times = times, unit = unit)
            # add the v, a, and w values to the dataframe
            mbm_res[row_idx, 1] <- RT[rt]
            mbm_res[row_idx, 2] <- V[v]
            mbm_res[row_idx, 3] <- A[a]
            mbm_res[row_idx, 4] <- W[w]
            mbm_res[row_idx, 5] <- SV[sv]
            # add the median microbenchmark results to the dataframe
            for (i in 1:nf) {
              mbm_res[row_idx, 5+i] <- median(mbm[mbm[,1] == fnames[i],2])
            }
            # iterate start value
            row_idx = row_idx + 1
          }
        }
      }
    }
  }
  return(mbm_res)
}




##### Parameter Space 1 Benchmarks #############################################
RT <- c(0.001, 0.1, 1, 2, 3, 4, 5, 10, 30)
A <- c(0.25, 0.5, 1, 2.5, 5)
V <- c(-5, -2, 0, 2, 5)
t0 <- 1e-4 # must be nonzero for RWiener
W <- seq(0.2, 0.8, by = 0.3)
SV <- seq(0, 1.5, by = 0.5)
err_tol <- 1e-6 # this is the setting from rtdists

bm_vec_1 <- rt_benchmark_vec(RT = RT, resp = "lower", V = V, A = A, t0 = t0,
                             W = W, SV = SV, err_tol = err_tol,
                             times = 1000, unit = "ns")

save(bm_vec_1, compress = "xz", compression_level = 9,
     file = paste0(save_dir, "mtl_bm_0-30.Rds"))

# uncomment the following line if loading pre-run benchmark data,
# will load into variable named "bm_vec_1"
# load(paste0(save_dir, "mtl_bm_0-30.Rds"))

# Plot Results
t_idx <- match("SV", colnames(bm_vec_1))
bm_vec_1[, -seq_len(t_idx)] <- bm_vec_1[, -seq_len(t_idx)]/1000 # convert to microseconds
mbm_vec_1 <- melt(bm_vec_1, measure.vars = -seq_len(t_idx),
                variable.name = "FuncName", value.name = "time")

Names_vec <- c("fb_SWSE_17_0", "fb_SWSE_14_0", "fb_SWSE_17_1", "fb_SWSE_14_1",
               "fb_SWSE_17_2", "fb_SWSE_14_2", "fb_SWSE_17_3", "fb_SWSE_14_3",
               "fb_SWSE_17_4", "fb_SWSE_14_4", "fb_SWSE_17_5", "fb_SWSE_14_5",
               "fb_SWSE_17_6", "fb_SWSE_14_6", "fb_SWSE_17_7", "fb_SWSE_14_7")
Color_vec <- c("#e000b4", "#ff99eb", "#e68a00", "#ffb366",
               "#006699", "#66ccff", "#9900cc", "#cc99ff",
               "#c2a500", "#d7db42", "#336600", "#33cc33",
               "#996633", "#ff9999", "#ff5050", "#990000")

ggplot(mbm_vec_1, aes(x = factor(FuncName, levels = Names_vec), y = time,
                    color = factor(FuncName, levels = Names_vec),
                    fill = factor(FuncName, levels = Names_vec))) +
  geom_violin(trim = TRUE, alpha = 0.5) +
  scale_color_manual(values = Color_vec, guide = FALSE) +
  scale_fill_manual(values = Color_vec, guide = FALSE) +
  geom_boxplot(width = 0.15, fill = "white", alpha = 0.5) +
  stat_summary(fun = mean, geom = "errorbar",
               aes(ymax = ..y.., ymin = ..y..),
               width = .35, linetype = "dashed") +
  labs(title = "Distribution of median benchmark times",
       subtitle = "Dashed lines represent mean benchmark times",
       x = "Method", y = "Time (microseconds)",
       color = "Method") +
  theme_bw() +
  theme(panel.border = element_blank(),
        plot.title = element_text(size = 23),
        plot.subtitle = element_text(size = 16),
        axis.text.x = element_text(size = 16, angle = 90,
                                   vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.position = "none")
ggsave(paste0(path, "mtl_bm_0-30.png"), width = 16, height = 9)

##### Parameter Space 2 Benchmorks #############################################
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
     file = paste0(save_dir, "mtl_bm_0-2.Rds"))

# uncomment the following line if loading pre-run benchmark data,
# will load into variable named "bm_vec_2"
# load(paste0(save_dir, "mtl_bm_0-2.Rds"))

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

ggplot(mbm_vec_2, aes(x = factor(FuncName, levels = Names_vec), y = time,
                    color = factor(FuncName, levels = Names_vec),
                    fill = factor(FuncName, levels = Names_vec))) +
  geom_violin(trim = TRUE, alpha = 0.5) +
  scale_color_manual(values = Color_vec, guide = FALSE) +
  scale_fill_manual(values = Color_vec, guide = FALSE) +
  geom_boxplot(width = 0.15, fill = "white", alpha = 0.5) +
  stat_summary(fun = mean, geom = "errorbar",
               aes(ymax = ..y.., ymin = ..y..),
               width = .35, linetype = "dashed") +
  labs(title = "Distribution of median benchmark times",
       subtitle = "Dashed lines represent mean benchmark times",
       x = "Method", y = "Time (microseconds)",
       color = "Method") +
  theme_bw() +
  theme(panel.border = element_blank(),
        plot.title = element_text(size = 23),
        plot.subtitle = element_text(size = 16),
        axis.text.x = element_text(size = 16, angle = 90,
                                   vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.position = "none")
ggsave(paste0(path, "mtl_bm_0-2.png"), width = 16, height = 9)
