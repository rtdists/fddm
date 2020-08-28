devtools::load_all(recompile = TRUE)
library("fddm")
library("microbenchmark")




rt_benchmark_vec <- function(RT, resp, V, A, t0 = 1e-4, W = 0.5, SV = 0.0,
                             err_tol = 1e-6, times = 100, unit = "ns") {

  fnames <- c("fc_SWSE_L1", "fc_SWSE_L2", "fc_SWSE_L3",
              "fc_SWSE_L4", "fd_SWSE_25", "fb_Gon_17")
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
            fc_SWSE_L1 = dfddm(rt = RT, response = resp, a = A[a],
                               v = V[v], t0 = t0, w = W[w],
                               log = FALSE, n_terms_small = "SWSE",
                               summation_small = "2017", scale = "both",
                               max_terms_large = 1, err_tol = err_tol),
            fc_SWSE_L2 = dfddm(rt = RT, response = resp, a = A[a],
                               v = V[v], t0 = t0, w = W[w],
                               log = FALSE, n_terms_small = "SWSE",
                               summation_small = "2017", scale = "both",
                               max_terms_large = 2, err_tol = err_tol),
            fc_SWSE_L3 = dfddm(rt = RT, response = resp, a = A[a],
                               v = V[v], t0 = t0, w = W[w],
                               log = FALSE, n_terms_small = "SWSE",
                               summation_small = "2017", scale = "both",
                               max_terms_large = 3, err_tol = err_tol),
            fc_SWSE_L4 = dfddm(rt = RT, response = resp, a = A[a],
                               v = V[v], t0 = t0, w = W[w],
                               log = FALSE, n_terms_small = "SWSE",
                               summation_small = "2017", scale = "both",
                               max_terms_large = 4, err_tol = err_tol),
            fd_SWSE_25 = dfddm(rt = RT, response = resp, a = A[a],
                               v = V[v], t0 = t0, w = W[w],
                               log = FALSE, n_terms_small = "SWSE",
                               summation_small = "2017", scale = "d",
                               max_terms_large = 5, err_tol = err_tol),
            fb_Gon_17 =  dfddm(rt = RT, response = resp, a = A[a],
                               v = V[v], t0 = t0, w = W[w],
                               log = FALSE, n_terms_small = "Gondan",
                               summation_small = "2017", scale = "both",
                               err_tol = err_tol),
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

  fnames <- c("fc_SWSE_L1", "fc_SWSE_L2", "fc_SWSE_L3",
              "fc_SWSE_L4", "fd_SWSE_25", "fb_Gon_17")
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
            res <- c(
              dfddm(rt = RT[rt], response = resp, a = A[a],
                    v = V[v], t0 = t0, w = W[w],
                    log = FALSE, n_terms_small = "SWSE",
                    summation_small = "2017", scale = "both",
                    max_terms_large = 1, err_tol = err_tol),
              dfddm(rt = RT[rt], response = resp, a = A[a],
                    v = V[v], t0 = t0, w = W[w],
                    log = FALSE, n_terms_small = "SWSE",
                    summation_small = "2017", scale = "both",
                    max_terms_large = 2, err_tol = err_tol),
              dfddm(rt = RT[rt], response = resp, a = A[a],
                    v = V[v], t0 = t0, w = W[w],
                    log = FALSE, n_terms_small = "SWSE",
                    summation_small = "2017", scale = "both",
                    max_terms_large = 3, err_tol = err_tol),
              dfddm(rt = RT[rt], response = resp, a = A[a],
                    v = V[v], t0 = t0, w = W[w],
                    log = FALSE, n_terms_small = "SWSE",
                    summation_small = "2017", scale = "both",
                    max_terms_large = 4, err_tol = err_tol),
              dfddm(rt = RT[rt], response = resp, a = A[a],
                    v = V[v], t0 = t0, w = W[w],
                    log = FALSE, n_terms_small = "SWSE",
                    summation_small = "2017", scale = "d",
                    max_terms_large = 5, err_tol = err_tol),
              dfddm(rt = RT[rt], response = resp, a = A[a],
                    v = V[v], t0 = t0, w = W[w],
                    log = FALSE, n_terms_small = "Gondan",
                    summation_small = "2017", scale = "both",
                    err_tol = err_tol)
            )
            mbm <- microbenchmark(
              fc_SWSE_L1 = dfddm(rt = RT[rt], response = resp, a = A[a],
                                 v = V[v], t0 = t0, w = W[w],
                                 log = FALSE, n_terms_small = "SWSE",
                                 summation_small = "2017", scale = "both",
                                 max_terms_large = 1, err_tol = err_tol),
              fc_SWSE_L2 = dfddm(rt = RT[rt], response = resp, a = A[a],
                                 v = V[v], t0 = t0, w = W[w],
                                 log = FALSE, n_terms_small = "SWSE",
                                 summation_small = "2017", scale = "both",
                                 max_terms_large = 2, err_tol = err_tol),
              fc_SWSE_L3 = dfddm(rt = RT[rt], response = resp, a = A[a],
                                 v = V[v], t0 = t0, w = W[w],
                                 log = FALSE, n_terms_small = "SWSE",
                                 summation_small = "2017", scale = "both",
                                 max_terms_large = 3, err_tol = err_tol),
              fc_SWSE_L4 = dfddm(rt = RT[rt], response = resp, a = A[a],
                                 v = V[v], t0 = t0, w = W[w],
                                 log = FALSE, n_terms_small = "SWSE",
                                 summation_small = "2017", scale = "both",
                                 max_terms_large = 4, err_tol = err_tol),
              fd_SWSE_25 = dfddm(rt = RT[rt], response = resp, a = A[a],
                                 v = V[v], t0 = t0, w = W[w],
                                 log = FALSE, n_terms_small = "SWSE",
                                 summation_small = "2017", scale = "d",
                                 max_terms_large = 5, err_tol = err_tol),
              fb_Gon_17 =  dfddm(rt = RT[rt], response = resp, a = A[a],
                                 v = V[v], t0 = t0, w = W[w],
                                 log = FALSE, n_terms_small = "Gondan",
                                 summation_small = "2017", scale = "both",
                                 err_tol = err_tol),
            times = times, unit = unit)
            # add the v, a, and w values to the dataframe
            mbm_res[row_idx, 1] <- RT[rt]
            mbm_res[row_idx, 2] <- V[v]
            mbm_res[row_idx, 3] <- A[a]
            mbm_res[row_idx, 4] <- W[w]
            mbm_res[row_idx, 5] <- SV[sv]
            # add the median microbenchmark results to the dataframe
            for (i in 1:nf) {
              tmp <- median(mbm[mbm[,1] == fnames[i],2])
              mbm_res[row_idx, 5+i] <- ifelse(res[i] < 0, -tmp, tmp)
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




RT <- c(0.001, 0.1, 1, 2, 3, 4, 5, 10, 30)
A <- c(0.25, 0.5, 1, 2.5, 5)
V <- c(-5, -2, 0, 2, 5)
t0 <- 1e-4 # must be nonzero for RWiener
W <- c(0.2, 0.5, 0.8)
SV <- c(0, 0.5, 1, 1.5)
err_tol <- 1e-6 # this is the setting from rtdists


bm_vec <- rt_benchmark_vec(RT = RT, resp = "lower", V = V, A = A, t0 = t0,
                           W = W, SV = SV, err_tol = err_tol,
                           times = 100, unit = "ns")

bm_ind <- rt_benchmark_ind(RT = RT, resp = "lower", V = V, A = A, t0 = t0,
                           W = W, SV = SV, err_tol = err_tol,
                           times = 100, unit = "ns")

save(bm_vec, compress = "xz", compression_level = 9,
     file = "zedonked/results/bm_vec.Rds")
save(bm_ind, compress = "xz", compression_level = 9,
     file = "zedonked/results/bm_ind.Rds")
load(file = "zedonked/results/bm_vec.Rds")
load(file = "zedonked/results/bm_ind.Rds")


library("reshape2")
library("ggplot2")

Names <- c("fc_SWSE_L1", "fc_SWSE_L2", "fc_SWSE_L3",
           "fc_SWSE_L4", "fd_SWSE_25", "fb_Gon_17")
Color <- c("#ff99fa", "#ff70f8", "#ff40f5",
           "#ff1ff3", "#39ff14", "#c2a500")

t_idx_vec <- match("SV", colnames(bm_vec))
bm_vec[, -seq_len(t_idx_vec)] <- bm_vec[, -seq_len(t_idx_vec)]/1000 # convert to microseconds
mbm_vec <- melt(bm_vec, measure.vars = -seq_len(t_idx_vec),
                variable.name = "FuncName", value.name = "time")
min_vec <- min(bm_vec[, -seq_len(t_idx_vec)])
max_vec <- max(bm_vec[, -seq_len(t_idx_vec)])

bm_ind$RTAA <- bm_ind$RT / bm_ind$A / bm_ind$A
bm_ind <- bm_ind[, c(1, 2, ncol(bm_ind), 3:(ncol(bm_ind)-1)) ]

t_idx_ind <- match("SV", colnames(bm_ind))
bm_ind[,-seq_len(t_idx_ind)] <- bm_ind[, -seq_len(t_idx_ind)]/1000 # convert to microseconds
mbm_ind <- melt(bm_ind, measure.vars = -seq_len(t_idx_ind),
                variable.name = "FuncName", value.name = "time")
mbm_ind[["Scale"]] <- ifelse(mbm_ind[["time"]] < 0, 0, 1)
mbm_ind[["time"]] <- abs(mbm_ind[["time"]])
min_ind <- min(mbm_ind[["time"]])
max_ind <- max(mbm_ind[["time"]])

mbm_ind_tail <- data.frame(matrix(ncol = ncol(mbm_ind), nrow = 0))
colnames(mbm_ind_tail) <- colnames(mbm_ind)
for (i in 1:length(Names)) {
  mbm_ind_tail <- rbind(mbm_ind_tail,
                        mbm_ind[mbm_ind[["FuncName"]] == Names[i] &
                                mbm_ind[["time"]] > median(mbm_ind[["time"]]) +
                                                    2*sd(mbm_ind[["time"]]), ])
}


# density (vec)
ggplot(mbm_vec, aes(x = FuncName, y = time,
                    color = factor(FuncName, levels = Names),
                    fill = factor(FuncName, levels = Names))) +
       geom_violin(trim = TRUE, alpha = 0.5) +
       scale_color_manual(values = Color) +
       scale_fill_manual(values = Color) +
       geom_boxplot(width = 0.15, fill = "white", alpha = 0.5) +
       stat_summary(fun = mean, geom = "errorbar",
                    aes(ymax = ..y.., ymin = ..y..),
                    width = .35, linetype = "dashed") +
       coord_cartesian(ylim = c(min_vec, max_vec)) +
       labs(title = "Distribution of median benchmark times (vec)",
            subtitle = "Dashed lines represent mean benchmark times",
            x = "Method", y = "Time (microseconds)") +
       theme_bw() +
       theme(panel.border = element_blank(),
             plot.title = element_text(size = 23),
             plot.subtitle = element_text(size = 16),
             axis.title.x = element_text(size = 18,
                                         margin = margin(10, 0, 0, 0, "pt")),
             axis.title.y = element_text(size = 18,
                                         margin = margin(0, 5, 0, 0, "pt")),
             axis.text.x = element_text(size = 16, angle = 90,
                                        vjust = 0.5, hjust = 1),
             axis.text.y = element_text(size = 16),
             legend.position = "none")
ggsave("zedonked/results/vectorized_rt_densities.png",
       width = 16, height = 9)

# density (ind)
ggplot(mbm_ind, aes(x = factor(FuncName, levels = Names), y = time,
                    color = factor(FuncName, levels = Names),
                    fill = factor(FuncName, levels = Names))) +
       geom_violin(trim = TRUE, alpha = 0.5) +
       scale_color_manual(values = Color, guide = FALSE) +
       scale_fill_manual(values = Color, guide = FALSE) +
       geom_boxplot(width = 0.15, fill = "white", alpha = 0.5,
                    outlier.shape = NA) +
       stat_summary(fun = mean, geom = "errorbar",
                    aes(ymax = ..y.., ymin = ..y..),
                    width = .35, linetype = "dashed") +
       geom_point(data = mbm_ind_tail,
                  aes(color = factor(FuncName, levels = Names),
                      size = factor(Scale, levels = c(0, 1)),
                      shape = factor(Scale, levels = c(0, 1)))) +
       scale_size_manual(values = c(3, 2), guide = FALSE) +
       scale_shape_manual(values = c(1, 4),
                          name = "Time Scale Density",
                          breaks = c(0, 1),
                          labels = c("Large Time", "Small Time")) +
       coord_cartesian(ylim = c(min_ind, max_ind)) +
       labs(title = "Distribution of median benchmark times (ind)",
            subtitle = paste("Dashed lines represent mean benchmark times;",
                             "Plotted points are in the tails",
                             sep = "\n"),
            x = "Method", y = "Time (microseconds)") +
       guides(shape = guide_legend(override.aes = list(size = c(3, 2)))) +
       theme_bw() +
       theme(panel.border = element_blank(),
             plot.title = element_text(size = 23),
             plot.subtitle = element_text(size = 16),
             axis.title.x = element_text(size = 18,
                                         margin = margin(10, 0, 0, 0, "pt")),
             axis.title.y = element_text(size = 18,
                                         margin = margin(0, 5, 0, 0, "pt")),
             axis.text.x = element_text(size = 16, angle = 90,
                                        vjust = 0.5, hjust = 1),
             axis.text.y = element_text(size = 16),
             legend.position = c(0, 1),
             legend.justification = c(0, 1),
             legend.box = "vertical",
             legend.direction = "vertical",
             legend.background = element_rect(fill = "transparent"),
             legend.title = element_text(size = 14),
             legend.text = element_text(size = 13))
ggsave("zedonked/results/individual_rt_densities.png",
       width = 16, height = 9)

# density (large time vs small time)
ggplot(mbm_ind[mbm_ind[["FuncName"]] == Names[2], ],
       aes(x = factor(Scale, levels = c(0, 1)), y = time,
           color = factor(Scale, levels = c(0, 1)),
           fill = factor(Scale, levels = c(0, 1)))) +
       geom_violin(trim = TRUE, alpha = 0.5) +
       scale_color_manual(values = c("#4287f5", "#F5B042"), guide = FALSE) +
       scale_fill_manual(values = c("#4287f5", "#F5B042"), guide = FALSE) +
       geom_boxplot(width = 0.15, fill = "white", alpha = 0.5,
                    outlier.shape = NA) +
       stat_summary(fun = mean, geom = "errorbar",
                    aes(ymax = ..y.., ymin = ..y..),
                    width = .35, linetype = "dashed") +
       coord_cartesian(ylim = c(min_ind, max_ind)) +
       scale_x_discrete(labels = c("0" = "Large Time", "1" = "Small Time")) +
       labs(title = "Distribution of median benchmark times (ind)",
            subtitle = paste(paste0("For each use of time scale in ",
                                    Names[2], ";"),
                             "Dashed lines represent mean benchmark times",
                             sep = "\n"),
            x = "Time Scale Method", y = "Time (microseconds)") +
       theme_bw() +
       theme(panel.border = element_blank(),
             plot.title = element_text(size = 23),
             plot.subtitle = element_text(size = 16),
             axis.title.x = element_text(size = 18,
                                         margin = margin(10, 0, 0, 0, "pt")),
             axis.title.y = element_text(size = 18,
                                         margin = margin(0, 5, 0, 0, "pt")),
             axis.text.x = element_text(size = 16),
             axis.text.y = element_text(size = 16),
             legend.position = "none")
ggsave("zedonked/results/large_vs_small_densities.png",
       width = 16, height = 9)

# mean with quantiles (ind)
ggplot(mbm_ind, aes(x = RTAA, y = time,
                    color = factor(FuncName, levels = Names),
                    fill = factor(FuncName, levels = Names))) +
  stat_summary(fun.min = min, fun.max = max,
               geom = "ribbon", color = NA, alpha = 0.1) +
  stat_summary(fun.min = function(z) { quantile(z, 0.1) },
               fun.max = function(z) { quantile(z, 0.9) },
               geom = "ribbon", color = NA, alpha = 0.2) +
  stat_summary(fun = mean, geom = "line") +
  geom_point(data = mbm_ind_tail,
             aes(x = RTAA, y = time,
                 color = factor(FuncName, levels = Names),
                 size = factor(Scale, levels = c(0, 1)),
                 shape = factor(Scale, levels = c(0, 1)))) +
  geom_vline(xintercept = 0.25, color = "grey60", alpha = 0.5) +
  scale_x_log10() +
  scale_color_manual(values = Color, guide = FALSE) +
  scale_fill_manual(values = Color, guide = FALSE) +
  scale_size_manual(values = c(2, 1.5), guide = FALSE) +
  scale_shape_manual(values = c(1, 2),
                     name = "Time Scale Density",
                     breaks = c(0, 1),
                     labels = c("Large Time", "Small Time"),
                     guide = guide_legend(direction = "horizontal",
                                          title.position = "top")) +
  labs(title = "Means of the median microbenchmark results",
       subtitle =
        paste("The shaded regions represent the 10% and 90% quantiles;",
              "The lighter shaded regions represent the min and max times;",
              "Plotted points are in the tails",
              sep = "\n"),
       x = bquote(frac(rt, a^2) ~ ", effective response time, " ~ log[10]),
       y = "Time (microseconds)") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        plot.title = element_text(size = 20),
        plot.subtitle = element_text(size = 16,
                                     margin = margin(5, 5, 40, 5, "pt")),
        axis.title.x = element_text(size = 18,
                                    margin = margin(5, 0, 0, 0, "pt")),
        axis.title.y = element_text(size = 18,
                                    margin = margin(0, 5, 0, 0, "pt")),
        axis.text.x = element_text(size = 16, angle = 90,
                                   vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 16),
        strip.text = element_text(size = 14),
        strip.background = element_rect(fill = "white"),
        legend.position = c(1, 1.1),
        legend.justification = c(1, 0),
        legend.background = element_rect(fill = "transparent"),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 13)) +
  facet_wrap(~ factor(FuncName, levels = Names))
ggsave("zedonked/results/eff_rt_mean_quantiles.png",
       width = 16, height = 9)
