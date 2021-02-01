# This file produces the results (and plots) that pertain to the raw algorithms,
# as shown in Appendix B.1 of the fddm paper: "Additional Results"

library("fddm")
library("rtdists")
library("RWiener")
source(system.file("extdata", "Gondan_et_al_density.R", package = "fddm", mustWork = TRUE))
library("microbenchmark")
library("reshape2")
library("ggplot2")
library("ggforce")
save_dir <- "paper_analysis/appendices/"



##### Benchmark Function #######################################################

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
RT <- seq(0.1, 2, by = 0.1)
A <- seq(0.5, 3.5, by = 0.5)
V <- c(-5, -2, 0, 2, 5)
t0 <- 1e-4 # must be nonzero for RWiener
W <- seq(0.3, 0.7, by = 0.1)
SV <- c(0, 1, 2, 3.5)
err_tol <- 1e-6 # this is the setting from rtdists

# NOTE: if already ran and saved benchmark data, can skip this
bm_ind <- rt_benchmark_ind(RT = RT, resp = "lower", V = V, A = A, t0 = t0,
                           W = W, SV = SV, err_tol = err_tol,
                           times = 100, unit = "ns")
save(bm_ind, compress = "xz", compression_level = 9,
     file = paste0(save_dir, "bm_ind_0-2.Rds"))



##### Plot Benchmark Timing Results ############################################
# uncomment the following line if loading pre-run benchmark data,
# will load into variable named "bm_ind"
# load("paper_analysis/section_5_benchmarking/raw_algorithms/bm_ind_0-2.Rds")

bm_ind[["RTAA"]] <- bm_ind[["RT"]] / bm_ind[["A"]] / bm_ind[["A"]]
bm_ind <- bm_ind[, c(1, ncol(bm_ind), 2:(ncol(bm_ind)-1)) ]
t_idx_ind <- match("SV", colnames(bm_ind))
bm_ind[, -seq_len(t_idx_ind)] <- bm_ind[, -seq_len(t_idx_ind)]/1000 # convert to microseconds
mbm_ind <- melt(bm_ind, measure.vars = -seq_len(t_idx_ind),
                variable.name = "FuncName", value.name = "time")

Names_ind <- c("fb_SWSE_17", "fb_SWSE_14", "fb_Gon_17", "fb_Gon_14",
               "fb_Nav_17", "fb_Nav_14", "fs_SWSE_17", "fs_SWSE_14",
               "fs_Gon_17", "fs_Gon_14", "fs_Nav_17", "fs_Nav_14",
               "fl_Nav_09", "RWiener", "Gondan", "rtdists")
Color_ind <- c("#e000b4", "#ff99eb", "#e68a00", "#ffb366",
               "#006699", "#66ccff", "#9900cc", "#cc99ff",
               "#c2a500", "#d7db42", "#336600", "#33cc33",
               "#996633", "#ff9999", "#ff5050", "#990000")

mi_ind <- min(bm_ind[, -seq_len(t_idx_ind)])
ma_ind <- max(bm_ind[, (t_idx_ind+1):(ncol(bm_ind)-4)])

ggplot(mbm_ind, aes(x = factor(FuncName, levels = Names_ind), y = time,
                    color = factor(FuncName, levels = Names_ind),
                    fill = factor(FuncName, levels = Names_ind))) +
       geom_violin(trim = TRUE, alpha = 0.5) +
       scale_color_manual(values = Color_ind, guide = FALSE) +
       scale_fill_manual(values = Color_ind, guide = FALSE) +
       geom_boxplot(width = 0.15, fill = "white", alpha = 0.5) +
       stat_summary(fun = mean, geom = "errorbar",
                    aes(ymax = ..y.., ymin = ..y..),
                    width = .35, linetype = "dashed") +
       facet_zoom(ylim = c(mi_ind, ma_ind)) +
       labs(title = "Distribution of median benchmark times",
            subtitle = paste("Response time data input as individual data points",
                             "Dashed lines represent mean benchmark times",
                             sep = ";\n"),
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
ggsave("paper_analysis/images/bm_ind_0-2.png", width = 16, height = 9)
