# Figures S3, S4, S5, S8, S9, S10

# This file produces the supplementary results (and plots) that pertain to all
# of the implementations discussed in the fddm paper on a predefined parameter
# space, as shown in
# Section 2.2 of the supplementary analysis document of the fddm paper:
# "Benchmarking All Implementations/Predefined Parameter Spaces"

library("fddm")
library("rtdists")
library("RWiener")
source(system.file("extdata", "Gondan_et_al_density.R", package = "fddm", mustWork = TRUE))
library("microbenchmark")
library("reshape2")
library("ggplot2")
library("ggforce")
save_dir <- "paper_analysis/supplemental_analysis/3_all_implementations/results/"
img_dir <- "paper_analysis/images/supplemental_analysis/3_all_implementations/"



##### Benchmark Functions ######################################################
rt_benchmark_vec <- function(RT, resp, V, A, t0 = 1e-4, W = 0.5, SV = 0.0,
                             err_tol = 1e-6, times = 100, unit = "ns") {

  fnames <- c("fs_SWSE_17", "fs_SWSE_14", "fb_SWSE_17", "fb_SWSE_14",
              "fs_Gon_17", "fs_Gon_14", "fb_Gon_17", "fb_Gon_14",
              "fs_Nav_17", "fs_Nav_14", "fb_Nav_17", "fb_Nav_14",
              "fl_Nav_09", "RWiener", "Gondan", "rtdists")
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
          fs_SWSE_17 = dfddm(rt = RT, response = resp, a = A[a],
                             v = V[v], t0 = t0, w = W[w],
                             log = FALSE, n_terms_small = "SWSE",
                             summation_small = "2017", scale = "small",
                             err_tol = err_tol),
          fs_SWSE_14 = dfddm(rt = RT, response = resp, a = A[a],
                             v = V[v], t0 = t0, w = W[w],
                             log = FALSE, n_terms_small = "SWSE",
                             summation_small = "2014", scale = "small",
                             err_tol = err_tol),
          fb_SWSE_17 = dfddm(rt = RT, response = resp, a = A[a],
                             v = V[v], t0 = t0, w = W[w],
                             log = FALSE, n_terms_small = "SWSE",
                             summation_small = "2017", scale = "both",
                             err_tol = err_tol),
          fb_SWSE_14 = dfddm(rt = RT, response = resp, a = A[a],
                             v = V[v], t0 = t0, w = W[w],
                             log = FALSE, n_terms_small = "SWSE",
                             summation_small = "2014", scale = "both",
                             err_tol = err_tol),
          fs_Gon_17 = dfddm(rt = RT, response = resp, a = A[a],
                            v = V[v], t0 = t0, w = W[w],
                            log = FALSE, n_terms_small = "Gondan",
                            summation_small = "2017", scale = "small",
                            err_tol = err_tol),
          fs_Gon_14 = dfddm(rt = RT, response = resp, a = A[a],
                            v = V[v], t0 = t0, w = W[w],
                            log = FALSE, n_terms_small = "Gondan",
                            summation_small = "2014", scale = "small",
                            err_tol = err_tol),
          fb_Gon_17 = dfddm(rt = RT, response = resp, a = A[a],
                            v = V[v], t0 = t0, w = W[w],
                            log = FALSE, n_terms_small = "Gondan",
                            summation_small = "2017", scale = "both",
                            err_tol = err_tol),
          fb_Gon_14 = dfddm(rt = RT, response = resp, a = A[a],
                            v = V[v], t0 = t0, w = W[w],
                            log = FALSE, n_terms_small = "Gondan",
                            summation_small = "2014", scale = "both",
                            err_tol = err_tol),
          fs_Nav_17 = dfddm(rt = RT, response = resp, a = A[a],
                            v = V[v], t0 = t0, w = W[w],
                            log = FALSE, n_terms_small = "Navarro",
                            summation_small = "2017", scale = "small",
                            err_tol = err_tol),
          fs_Nav_14 = dfddm(rt = RT, response = resp, a = A[a],
                            v = V[v], t0 = t0, w = W[w],
                            log = FALSE, n_terms_small = "Navarro",
                            summation_small = "2014", scale = "small",
                            err_tol = err_tol),
          fb_Nav_17 = dfddm(rt = RT, response = resp, a = A[a],
                            v = V[v], t0 = t0, w = W[w],
                            log = FALSE, n_terms_small = "Navarro",
                            summation_small = "2017", scale = "both",
                            err_tol = err_tol),
          fb_Nav_14 = dfddm(rt = RT, response = resp, a = A[a],
                            v = V[v], t0 = t0, w = W[w],
                            log = FALSE, n_terms_small = "Navarro",
                            summation_small = "2014", scale = "both",
                            err_tol = err_tol),
          fl_Nav_09 = dfddm(rt = RT, response = resp, a = A[a],
                            v = V[v], t0 = t0, w = W[w],
                            log = FALSE, n_terms_small = "Navarro",
                            scale = "large", err_tol = err_tol),
          RWiener = dwiener(RT, resp = resp, alpha = A[a],
                            delta = V[v], tau = t0, beta = W[w],
                            give_log = FALSE),
          Gondan = fs(t = RT-t0, a = A[a], v = V[v],
                      w = W[w], eps = err_tol), # only "lower" resp
          rtdists = ddiffusion(RT, resp, a = A[a], v = V[v],
                               t0 = t0, z = W[w]*A[a]),
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
  fnames <- c("fs_SWSE_17", "fs_SWSE_14", "fb_SWSE_17", "fb_SWSE_14",
              "fs_Gon_17", "fs_Gon_14", "fb_Gon_17", "fb_Gon_14",
              "fs_Nav_17", "fs_Nav_14", "fb_Nav_17", "fb_Nav_14",
              "fl_Nav_09", "RWiener", "Gondan", "rtdists")
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
            fs_SWSE_17 = dfddm(rt = RT[rt], response = resp, a = A[a],
                              v = V[v], t0 = t0, w = W[w],
                              log = FALSE, n_terms_small = "SWSE",
                              summation_small = "2017", scale = "small",
                              err_tol = err_tol),
            fs_SWSE_14 = dfddm(rt = RT[rt], response = resp, a = A[a],
                              v = V[v], t0 = t0, w = W[w],
                              log = FALSE, n_terms_small = "SWSE",
                              summation_small = "2014", scale = "small",
                              err_tol = err_tol),
            fb_SWSE_17 = dfddm(rt = RT[rt], response = resp, a = A[a],
                              v = V[v], t0 = t0, w = W[w],
                              log = FALSE, n_terms_small = "SWSE",
                              summation_small = "2017", scale = "both",
                              err_tol = err_tol),
            fb_SWSE_14 = dfddm(rt = RT[rt], response = resp, a = A[a],
                              v = V[v], t0 = t0, w = W[w],
                              log = FALSE, n_terms_small = "SWSE",
                              summation_small = "2014", scale = "both",
                              err_tol = err_tol),
            fs_Gon_17 = dfddm(rt = RT[rt], response = resp, a = A[a],
                              v = V[v], t0 = t0, w = W[w],
                              log = FALSE, n_terms_small = "Gondan",
                              summation_small = "2017", scale = "small",
                              err_tol = err_tol),
            fs_Gon_14 = dfddm(rt = RT[rt], response = resp, a = A[a],
                              v = V[v], t0 = t0, w = W[w],
                              log = FALSE, n_terms_small = "Gondan",
                              summation_small = "2014", scale = "small",
                              err_tol = err_tol),
            fb_Gon_17 = dfddm(rt = RT[rt], response = resp, a = A[a],
                              v = V[v], t0 = t0, w = W[w],
                              log = FALSE, n_terms_small = "Gondan",
                              summation_small = "2017", scale = "both",
                              err_tol = err_tol),
            fb_Gon_14 = dfddm(rt = RT[rt], response = resp, a = A[a],
                              v = V[v], t0 = t0, w = W[w],
                              log = FALSE, n_terms_small = "Gondan",
                              summation_small = "2014", scale = "both",
                              err_tol = err_tol),
            fs_Nav_17 = dfddm(rt = RT[rt], response = resp, a = A[a],
                              v = V[v], t0 = t0, w = W[w],
                              log = FALSE, n_terms_small = "Navarro",
                              summation_small = "2017", scale = "small",
                              err_tol = err_tol),
            fs_Nav_14 = dfddm(rt = RT[rt], response = resp, a = A[a],
                              v = V[v], t0 = t0, w = W[w],
                              log = FALSE, n_terms_small = "Navarro",
                              summation_small = "2014", scale = "small",
                              err_tol = err_tol),
            fb_Nav_17 = dfddm(rt = RT[rt], response = resp, a = A[a],
                              v = V[v], t0 = t0, w = W[w],
                              log = FALSE, n_terms_small = "Navarro",
                              summation_small = "2017", scale = "both",
                              err_tol = err_tol),
            fb_Nav_14 = dfddm(rt = RT[rt], response = resp, a = A[a],
                              v = V[v], t0 = t0, w = W[w],
                              log = FALSE, n_terms_small = "Navarro",
                              summation_small = "2014", scale = "both",
                              err_tol = err_tol),
            fl_Nav_09 = dfddm(rt = RT[rt], response = resp, a = A[a],
                              v = V[v], t0 = t0, w = W[w],
                              log = FALSE, n_terms_small = "Navarro",
                              scale = "large", err_tol = err_tol),
            RWiener = dwiener(RT[rt], resp = resp, alpha = A[a],
                              delta = V[v], tau = t0, beta = W[w],
                              give_log = FALSE),
            Gondan = fs(t = RT[rt]-t0, a = A[a], v = V[v],
                        w = W[w], eps = err_tol), # only "lower" resp
            rtdists = ddiffusion(RT[rt], resp, a = A[a], v = V[v],
                                 t0 = t0, z = W[w]*A[a]),
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



##### Run and Save the Benchmark Results #######################################
# Parameter Space RT range: 0-2 seconds
RT <- seq(0.1, 2, by = 0.1)
A <- seq(0.5, 3.5, by = 0.5)
V <- c(-5, -2, 0, 2, 5)
t0 <- 1e-4 # must be nonzero for RWiener
W <- seq(0.3, 0.7, by = 0.1)
SV <- c(0, 1, 2, 3.5)
err_tol <- 1e-6 # this is the setting from rtdists

bm_ind <- rt_benchmark_ind(RT = RT, resp = "lower", V = V, A = A, t0 = t0,
                           W = W, SV = SV, err_tol = err_tol,
                           times = 1000, unit = "ns")
save(bm_ind, compress = "xz", compression_level = 9,
     file = paste0(save_dir, "bm_ind_0-2.Rds"))


# Parameter Space RT range: 0-30 seconds
RT <- c(0.001, 0.1, 1, 2, 3, 4, 5, 10, 30)
A <- c(0.25, 0.5, 1, 2.5, 5)
V <- c(-5, -2, 0, 2, 5)
t0 <- 1e-4 # must be nonzero for RWiener
W <- c(0.2, 0.5, 0.8)
SV <- c(0, 0.5, 1, 1.5)
err_tol <- 1e-6 # this is the setting from rtdists
bm_vec <- rt_benchmark_vec(RT = RT, resp = "lower", V = V, A = A, t0 = t0,
                           W = W, SV = SV, err_tol = err_tol,
                           times = 10000, unit = "ns")
save(bm_vec, compress = "xz", compression_level = 9,
     file = paste0(save_dir, "bm_vec_0-30.Rds"))

bm_ind <- rt_benchmark_ind(RT = RT, resp = "lower", V = V, A = A, t0 = t0,
                           W = W, SV = SV, err_tol = err_tol,
                           times = 1000, unit = "ns")
save(bm_ind, compress = "xz", compression_level = 9,
     file = paste0(save_dir, "bm_ind_0-30.Rds"))



##### Plot Benchmark Timing Results ############################################

##### Vectorized Inputs (RT)
# uncomment the following line if loading pre-run benchmark data,
# will load into variable named "bm_vec"
# load(paste0(save_dir, "bm_vec_0-30.Rds"))

t_idx <- match("SV", colnames(bm_vec))
bm_vec[, -seq_len(t_idx)] <- bm_vec[, -seq_len(t_idx)]/1000 # convert to microseconds
mbm_vec <- melt(bm_vec, measure.vars = -seq_len(t_idx),
                variable.name = "FuncName", value.name = "time")

Names_vec <- c("fb_SWSE_17", "fb_SWSE_14", "fb_Gon_17", "fb_Gon_14",
               "fb_Nav_17", "fb_Nav_14", "fs_SWSE_17", "fs_SWSE_14",
               "fs_Gon_17", "fs_Gon_14", "fs_Nav_17", "fs_Nav_14",
               "fl_Nav_09", "RWiener", "Gondan", "rtdists")
Color_vec <- c("#92c639", "#d3e8b0", "#b3724d", "#e0c7b8",
               "#4da7b3", "#b8dce0", "#5cc639", "#bee8b0",
               "#b34d4d", "#e0b8b8", "#4d80b3", "#b8cce0",
               "#dcdca3", "#deccba", "#c5a687", "#ac8053")
Outline_vec <- c("#92c639", "#92c639", "#b3724d", "#b3724d",
               "#4da7b3", "#4da7b3", "#5cc639", "#5cc639",
               "#b34d4d", "#b34d4d", "#4d80b3", "#4d80b3",
               "#dcdca3", "#deccba", "#c5a687", "#ac8053")

mi <- min(bm_vec[, -seq_len(t_idx)])
ma <- max(bm_vec[, (t_idx+1):(ncol(bm_vec)-4)])

# Figure S4
fig_s4 <- ggplot(mbm_vec, aes(x = factor(FuncName, levels = Names_vec), y = time,
                              color = factor(FuncName, levels = Names_vec),
                              fill = factor(FuncName, levels = Names_vec))) +
  geom_violin(trim = TRUE, alpha = 0.5) +
  scale_color_manual(values = Outline_vec, guide = FALSE) +
  scale_fill_manual(values = Color_vec, guide = FALSE) +
  geom_boxplot(width = 0.15, fill = "white", alpha = 0.5) +
  stat_summary(fun = mean, geom = "errorbar",
               aes(ymax = ..y.., ymin = ..y..),
               width = .35, linetype = "dashed") +
  scale_x_discrete(labels = c(
    bquote(f[c] ~ SWSE[17]), bquote(f[c] ~ SWSE[14]),
    bquote(f[c] ~ Gon[17]), bquote(f[c] ~ Gon[14]),
    bquote(f[c] ~ Nav[17]), bquote(f[c] ~ Nav[14]),
    bquote(f[s] ~ SWSE[17]), bquote(f[s] ~ SWSE[14]),
    bquote(f[s] ~ Gon[17]), bquote(f[s] ~ Gon[14]),
    bquote(f[s] ~ Nav[17]), bquote(f[s] ~ Nav[14]),
    bquote(f[l] ~ "Nav"), "RWiener",
    "Gondan", "rtdists")) +
  facet_zoom(ylim = c(mi, ma)) +
  labs(x = "Implementation", y = "Time (microseconds)") +
  theme_bw() +
  theme(panel.border = element_blank(),
        axis.text.x = element_text(size = 16, angle = 90,
                                   vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 20,
                                    margin = margin(10, 0, 0, 0)),
        axis.title.y = element_text(size = 20,
                                    margin = margin(0, 10, 0, 0)),
        legend.position = "none")
ggsave(paste0(img_dir, "bm_vec_0-30.png"),
       plot = fig_s4, width = 16, height = 9)


##### Non-Vectorized Inputs (RT)
# uncomment the following line if loading pre-run benchmark data,
# will load into variable named "bm_ind"
# load(paste0(save_dir, "bm_ind_0-2.Rds"))

t_idx <- match("SV", colnames(bm_ind))
bm_ind[, -seq_len(t_idx)] <- bm_ind[, -seq_len(t_idx)]/1000 # convert to microseconds
mbm_ind <- melt(bm_ind, measure.vars = -seq_len(t_idx),
                variable.name = "FuncName", value.name = "time")

Names_ind <- c("fb_SWSE_17", "fb_SWSE_14", "fb_Gon_17", "fb_Gon_14",
               "fb_Nav_17", "fb_Nav_14", "fs_SWSE_17", "fs_SWSE_14",
               "fs_Gon_17", "fs_Gon_14", "fs_Nav_17", "fs_Nav_14",
               "fl_Nav_09", "RWiener", "Gondan", "rtdists")
Color_ind <- c("#92c639", "#d3e8b0", "#b3724d", "#e0c7b8",
               "#4da7b3", "#b8dce0", "#5cc639", "#bee8b0",
               "#b34d4d", "#e0b8b8", "#4d80b3", "#b8cce0",
               "#dcdca3", "#deccba", "#c5a687", "#ac8053")
Outline_ind <- c("#92c639", "#92c639", "#b3724d", "#b3724d",
               "#4da7b3", "#4da7b3", "#5cc639", "#5cc639",
               "#b34d4d", "#b34d4d", "#4d80b3", "#4d80b3",
               "#dcdca3", "#deccba", "#c5a687", "#ac8053")


mi <- min(bm_ind[, -seq_len(t_idx)])
ma <- max(bm_ind[, (t_idx+1):(ncol(bm_ind)-4)])

# Figure S3
fig_s3 <- ggplot(mbm_ind, aes(x = factor(FuncName, levels = Names_ind), y = time,
                              color = factor(FuncName, levels = Names_ind),
                              fill = factor(FuncName, levels = Names_ind))) +
  geom_violin(trim = TRUE, alpha = 0.5) +
  scale_color_manual(values = Outline_ind, guide = FALSE) +
  scale_fill_manual(values = Color_ind, guide = FALSE) +
  geom_boxplot(width = 0.15, fill = "white", alpha = 0.5) +
  stat_summary(fun = mean, geom = "errorbar",
               aes(ymax = ..y.., ymin = ..y..),
               width = .35, linetype = "dashed") +
  scale_x_discrete(labels = c(
    bquote(f[c] ~ SWSE[17]), bquote(f[c] ~ SWSE[14]),
    bquote(f[c] ~ Gon[17]), bquote(f[c] ~ Gon[14]),
    bquote(f[c] ~ Nav[17]), bquote(f[c] ~ Nav[14]),
    bquote(f[s] ~ SWSE[17]), bquote(f[s] ~ SWSE[14]),
    bquote(f[s] ~ Gon[17]), bquote(f[s] ~ Gon[14]),
    bquote(f[s] ~ Nav[17]), bquote(f[s] ~ Nav[14]),
    bquote(f[l] ~ "Nav"), "RWiener",
    "Gondan", "rtdists")) +
  facet_zoom(ylim = c(mi, ma)) +
  labs(x = "Implementation", y = "Time (microseconds)") +
  theme_bw() +
  theme(panel.border = element_blank(),
        axis.text.x = element_text(size = 16, angle = 90,
                                   vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 20,
                                    margin = margin(10, 0, 0, 0)),
        axis.title.y = element_text(size = 20,
                                    margin = margin(0, 10, 0, 0)),
        legend.position = "none")
ggsave(paste0(img_dir, "bm_ind_0-2.png"),
       plot = fig_s3, width = 16, height = 9)



# uncomment the following line if loading pre-run benchmark data,
# will load into variable named "bm_ind"
# load(paste0(save_dir, "bm_ind_0-30.Rds"))
bm_ind[["RTAA"]] <- bm_ind[["RT"]] / bm_ind[["A"]] / bm_ind[["A"]]
bm_ind <- bm_ind[, c(1, 2, ncol(bm_ind), 3:(ncol(bm_ind)-1)) ]

t_idx <- match("SV", colnames(bm_ind))
bm_ind[,-seq_len(t_idx)] <- bm_ind[, -seq_len(t_idx)]/1000 # convert to microseconds
mbm_ind <- melt(bm_ind, measure.vars = -seq_len(t_idx),
                variable.name = "FuncName", value.name = "time")


# Figure S5
fig_s5 <- ggplot(mbm_ind, aes(x = factor(FuncName, levels = Names_ind), y = time,
                              color = factor(FuncName, levels = Names_ind),
                              fill = factor(FuncName, levels = Names_ind))) +
  geom_violin(trim = TRUE, alpha = 0.5) +
  scale_color_manual(values = Outline_ind, guide = FALSE) +
  scale_fill_manual(values = Color_ind, guide = FALSE) +
  geom_boxplot(width = 0.15, fill = "white", alpha = 0.5) +
  stat_summary(fun = mean, geom = "errorbar",
               aes(ymax = ..y.., ymin = ..y..),
               width = .35, linetype = "dashed") +
  scale_x_discrete(labels = c(
    bquote(f[c] ~ SWSE[17]), bquote(f[c] ~ SWSE[14]),
    bquote(f[c] ~ Gon[17]), bquote(f[c] ~ Gon[14]),
    bquote(f[c] ~ Nav[17]), bquote(f[c] ~ Nav[14]),
    bquote(f[s] ~ SWSE[17]), bquote(f[s] ~ SWSE[14]),
    bquote(f[s] ~ Gon[17]), bquote(f[s] ~ Gon[14]),
    bquote(f[s] ~ Nav[17]), bquote(f[s] ~ Nav[14]),
    bquote(f[l] ~ "Nav"), "RWiener",
    "Gondan", "rtdists")) +
  facet_zoom(ylim = c(mi, ma)) +
  labs(x = "Implementation", y = "Time (microseconds)") +
  theme_bw() +
  theme(panel.border = element_blank(),
        axis.text.x = element_text(size = 16, angle = 90,
                                   vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 20,
                                    margin = margin(10, 0, 0, 0)),
        axis.title.y = element_text(size = 20,
                                    margin = margin(0, 10, 0, 0)),
        legend.position = "none")
ggsave(paste0(img_dir, "bm_ind_0-30.png"),
       plot = fig_s5, width = 16, height = 9)



### Mean and Quantile Plots
Names_meq <- c("fb_SWSE_17", "fs_SWSE_17", "fl_Nav_09",
               "RWiener", "Gondan", "rtdists")
Color_meq <- c("#92c639", "#5cc639", "#dcdca3",
               "#deccba", "#c5a687", "#ac8053")
mbm_meq <- subset(mbm_ind, FuncName %in% Names_meq)

my_labeller <- as_labeller(c(fb_SWSE_17 = "f[c] ~ SWSE[17]",
                             fs_SWSE_17 = "f[s] ~ SWSE[17]",
                             fl_Nav_09 = "f[l] ~ Nav",
                             RWiener = "RWiener",
                             Gondan = "Gondan",
                             rtdists = "rtdists"),
                           default = label_parsed)

### W (Relative Starting Point)
# Figure S8
fig_s8 <- ggplot(mbm_meq, aes(x = W, y = time,
                              color = factor(FuncName, levels = Names_meq),
                              fill = factor(FuncName, levels = Names_meq))) +
  stat_summary(fun.min = min, fun.max = max,
               geom = "ribbon", color = NA, alpha = 0.1) +
  stat_summary(fun.min = function(z) { quantile(z, 0.1) },
               fun.max = function(z) { quantile(z, 0.9) },
               geom = "ribbon", color = NA, alpha = 0.2) +
  stat_summary(fun = mean, geom = "line") +
  scale_color_manual(values = Color_meq) +
  scale_fill_manual(values = Color_meq) +
  labs(subtitle = paste(
         "The shaded regions represent the 10% and 90% quantiles",
         "The lighter shaded regions represent the min and max times",
         sep = ";\n"),
       x = "w, relative starting point (a priori bias)",
       y = "Time (microseconds)") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        plot.subtitle = element_text(size = 16,
                                     margin = margin(0, 0, 15, 0)),
        axis.text.x = element_text(size = 16, angle = 90,
                                   vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 20,
                                    margin = margin(10, 0, 0, 0)),
        axis.title.y = element_text(size = 20,
                                    margin = margin(0, 10, 0, 0)),
        strip.text = element_text(size = 16),
        strip.background = element_rect(fill = "white"),
        legend.position = "none") +
  facet_wrap(~ factor(FuncName, levels = Names_meq), scales = "free_y",
             labeller = my_labeller)
ggsave(paste0(img_dir, "meq_w_0-30.png"),
       plot = fig_s7, width = 16, height = 9)


### V (Drift Rate)
# Figure S9
fig_s9 <- ggplot(mbm_meq, aes(x = V, y = time,
                              color = factor(FuncName, levels = Names_meq),
                              fill = factor(FuncName, levels = Names_meq))) +
  stat_summary(fun.min = min, fun.max = max,
               geom = "ribbon", color = NA, alpha = 0.1) +
  stat_summary(fun.min = function(z) { quantile(z, 0.1) },
               fun.max = function(z) { quantile(z, 0.9) },
               geom = "ribbon", color = NA, alpha = 0.2) +
  stat_summary(fun = mean, geom = "line") +
  scale_color_manual(values = Color_meq) +
  scale_fill_manual(values = Color_meq) +
  labs(subtitle = paste(
         "The shaded regions represent the 10% and 90% quantiles",
         "The lighter shaded regions represent the min and max times",
         sep = ";\n"),
       x = "v, drift rate",
       y = "Time (microseconds)") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        plot.subtitle = element_text(size = 16,
                                     margin = margin(0, 0, 15, 0)),
        axis.text.x = element_text(size = 16, angle = 90,
                                   vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 20,
                                    margin = margin(10, 0, 0, 0)),
        axis.title.y = element_text(size = 20,
                                    margin = margin(0, 10, 0, 0)),
        strip.text = element_text(size = 16),
        strip.background = element_rect(fill = "white"),
        legend.position = "none") +
  facet_wrap(~ factor(FuncName, levels = Names_meq), scales = "free_y",
             labeller = my_labeller)
ggsave(paste0(img_dir, "meq_v_0-30.png"),
       plot = fig_s9, width = 16, height = 9)


### SV (Inter-Trial Variability in the Drift Rate)
# Figure S10
fig_s10 <- ggplot(mbm_meq, aes(x = SV, y = time,
                              color = factor(FuncName, levels = Names_meq),
                              fill = factor(FuncName, levels = Names_meq))) +
  stat_summary(fun.min = min, fun.max = max,
               geom = "ribbon", color = NA, alpha = 0.1) +
  stat_summary(fun.min = function(z) { quantile(z, 0.1) },
               fun.max = function(z) { quantile(z, 0.9) },
               geom = "ribbon", color = NA, alpha = 0.2) +
  stat_summary(fun = mean, geom = "line") +
  scale_color_manual(values = Color_meq) +
  scale_fill_manual(values = Color_meq) +
  labs(subtitle = paste(
         "The shaded regions represent the 10% and 90% quantiles",
         "The lighter shaded regions represent the min and max times",
         sep = ";\n"),
       x = "sv, inter-trial variability in the drift rate",
       y = "Time (microseconds)") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        plot.subtitle = element_text(size = 16,
                                     margin = margin(0, 0, 15, 0)),
        axis.text.x = element_text(size = 16, angle = 90,
                                   vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 20,
                                    margin = margin(10, 0, 0, 0)),
        axis.title.y = element_text(size = 20,
                                    margin = margin(0, 10, 0, 0)),
        strip.text = element_text(size = 16),
        strip.background = element_rect(fill = "white"),
        legend.position = "none") +
  facet_wrap(~ factor(FuncName, levels = Names_meq), scales = "free_y",
             labeller = my_labeller)
ggsave(paste0(img_dir, "meq_sv_0-30.png"),
       plot = fig_s10, width = 16, height = 9)
