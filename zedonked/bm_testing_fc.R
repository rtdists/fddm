devtools::load_all(recompile = TRUE)
library("fddm")
library("microbenchmark")




rt_benchmark_vec <- function(RT, resp, V, A, t0 = 1e-4, W = 0.5, SV = 0.0,
                             err_tol = 1e-6, times = 100, unit = "ns") {

  fnames <- c("L1", "L2", "L3", "L4", "L5")
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
            L1 = dfddm(rt = RT, response = resp, a = A[a],
                      v = V[v], t0 = t0, w = W[w],
                      log = FALSE, n_terms_small = "SWSE",
                      summation_small = "2017", scale = "both",
                      max_terms_large = 1, err_tol = err_tol),
            L2 = dfddm(rt = RT, response = resp, a = A[a],
                      v = V[v], t0 = t0, w = W[w],
                      log = FALSE, n_terms_small = "SWSE",
                      summation_small = "2017", scale = "both",
                      max_terms_large = 2, err_tol = err_tol),
            L3 = dfddm(rt = RT, response = resp, a = A[a],
                      v = V[v], t0 = t0, w = W[w],
                      log = FALSE, n_terms_small = "SWSE",
                      summation_small = "2017", scale = "both",
                      max_terms_large = 3, err_tol = err_tol),
            L4 = dfddm(rt = RT, response = resp, a = A[a],
                      v = V[v], t0 = t0, w = W[w],
                      log = FALSE, n_terms_small = "SWSE",
                      summation_small = "2017", scale = "both",
                      max_terms_large = 4, err_tol = err_tol),
            L5 = dfddm(rt = RT, response = resp, a = A[a],
                      v = V[v], t0 = t0, w = W[w],
                      log = FALSE, n_terms_small = "SWSE",
                      summation_small = "2017", scale = "both",
                      max_terms_large = 5, err_tol = err_tol),
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

  fnames <- c("L1", "L2", "L3", "L4", "L5")
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
              L1 = dfddm(rt = RT[rt], response = resp, a = A[a],
                        v = V[v], t0 = t0, w = W[w],
                        log = FALSE, n_terms_small = "SWSE",
                        summation_small = "2017", scale = "both",
                        max_terms_large = 1, err_tol = err_tol),
              L2 = dfddm(rt = RT[rt], response = resp, a = A[a],
                        v = V[v], t0 = t0, w = W[w],
                        log = FALSE, n_terms_small = "SWSE",
                        summation_small = "2017", scale = "both",
                        max_terms_large = 2, err_tol = err_tol),
              L3 = dfddm(rt = RT[rt], response = resp, a = A[a],
                        v = V[v], t0 = t0, w = W[w],
                        log = FALSE, n_terms_small = "SWSE",
                        summation_small = "2017", scale = "both",
                        max_terms_large = 3, err_tol = err_tol),
              L4 = dfddm(rt = RT[rt], response = resp, a = A[a],
                        v = V[v], t0 = t0, w = W[w],
                        log = FALSE, n_terms_small = "SWSE",
                        summation_small = "2017", scale = "both",
                        max_terms_large = 4, err_tol = err_tol),
              L5 = dfddm(rt = RT[rt], response = resp, a = A[a],
                        v = V[v], t0 = t0, w = W[w],
                        log = FALSE, n_terms_small = "SWSE",
                        summation_small = "2017", scale = "both",
                        max_terms_large = 5, err_tol = err_tol),
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




RT <- c(0.001, 0.1, 10, 30)
A <- c(0.5, 1)
V <- c(-2, 0, 2)
t0 <- 1e-4 # must be nonzero for RWiener
W <- c(0.2, 0.5, 0.8)
SV <- c(0, 1.5)
err_tol <- 1e-6 # this is the setting from rtdists


bm_vec <- rt_benchmark_vec(RT = RT, resp = "lower", V = V, A = A, t0 = t0,
                           W = W, SV = SV, err_tol = err_tol,
                           times = 100, unit = "ns")

bm_ind <- rt_benchmark_ind(RT = RT, resp = "lower", V = V, A = A, t0 = t0,
                           W = W, SV = SV, err_tol = err_tol,
                           times = 100, unit = "ns")




library("reshape2")
library("ggplot2")

t_idx <- match("SV", colnames(bm_vec))
bm_vec[, -seq_len(t_idx)] <- bm_vec[, -seq_len(t_idx)]/1000 # convert to microseconds
mbm_vec <- melt(bm_vec, measure.vars = -seq_len(t_idx),
                variable.name = "FuncName", value.name = "time")

Names_vec <- c("L1", "L2", "L3", "L4", "L5")
Color_vec <- c("#ff99fa", "#c2a500", "#e68a00",
               "#9900cc", "#006699")

mi <- min(bm_vec[, -seq_len(t_idx)])
ma <- max(bm_vec[, -seq_len(t_idx)])

ggplot(mbm_vec, aes(x = FuncName, y = time,
                    color = factor(FuncName, levels = Names_vec),
                    fill = factor(FuncName, levels = Names_vec))) +
       geom_violin(trim = TRUE, alpha = 0.5) +
       scale_color_manual(values = Color_vec) +
       scale_fill_manual(values = Color_vec) +
       geom_boxplot(width = 0.15, fill = "white", alpha = 0.5) +
       stat_summary(fun = mean, geom = "errorbar",
                    aes(ymax = ..y.., ymin = ..y..),
                    width = .35, linetype = "dashed") +
       coord_cartesian(ylim = c(mi, ma)) +
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





bm_ind$RTAA <- bm_ind$RT / bm_ind$A / bm_ind$A
bm_ind <- bm_ind[, c(1, 2, ncol(bm_ind), 3:(ncol(bm_ind)-1)) ]

t_idx <- match("SV", colnames(bm_ind))
bm_ind[,-seq_len(t_idx)] <- bm_ind[, -seq_len(t_idx)]/1000 # convert to microseconds
mbm_ind <- melt(bm_ind, measure.vars = -seq_len(t_idx),
                variable.name = "FuncName", value.name = "time")

Names_meq <- c("L1", "L2", "L3", "L4", "L5")
Color_meq <- c("#ff99fa", "#c2a500", "#e68a00",
               "#9900cc", "#006699")
mbm_meq <- subset(mbm_ind, FuncName %in% Names_meq)

ggplot(mbm_meq, aes(x = RTAA, y = time,
                    color = factor(FuncName, levels = Names_meq),
                    fill = factor(FuncName, levels = Names_meq))) +
  stat_summary(fun.min = min, fun.max = max,
               geom = "ribbon", color = NA, alpha = 0.1) +
  stat_summary(fun.min = function(z) { quantile(z, 0.1) },
               fun.max = function(z) { quantile(z, 0.9) },
               geom = "ribbon", color = NA, alpha = 0.2) +
  stat_summary(fun = mean, geom = "line") +
  scale_x_log10() +
  scale_color_manual(values = Color_meq) +
  scale_fill_manual(values = Color_meq) +
  labs(title = "Means of the median microbenchmark results",
       subtitle = "The shaded regions represent the 10% and 90% quantiles;\nThe lighter shaded regions represent the min and max times",
       x = bquote(frac(rt, a^2) ~ ", effective response time, " ~ log[10]),
       y = "Time (microseconds)") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        plot.title = element_text(size = 20),
        plot.subtitle = element_text(size = 16,
                                     margin = margin(5,5,15,5, "pt")),
        axis.text.x = element_text(size = 16, angle = 90,
                                   vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        strip.text = element_text(size = 14),
        strip.background = element_rect(fill = "white"),
        legend.position = "none") +
  facet_wrap(~ factor(FuncName, levels = Names_meq), scales = "free_y")
