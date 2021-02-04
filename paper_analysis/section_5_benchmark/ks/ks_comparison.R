library("Rcpp")
library("microbenchmark")
library("reshape2")
library("ggplot2")
save_dir <- "paper_analysis/section_5_benchmark/ks/"

# Load the ks functions (from C++, using Rcpp)
# the functions are:
  # ks_Gon(t, w, eps)
  # ks_Nav(t, w, eps) (NOTE: ks_Nav does not use the input `w`)
sourceCpp(paste0(save_dir, "ks.cpp"))



##### Benchmark Function #######################################################
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



##### Run and Save the Benchmark Results #######################################
ks_bm <- ks_benchmark(t_range = c(0, 2), dt = 0.02, def = 1000,
                      w = 0.5, eps = 1e-6, times = 1000)
save(ks_bm, compress = "xz", compression_level = 9,
     file = paste0(save_dir, "ks_bm.Rds"))



##### Plot Benchmark Timing Results ############################################
# uncomment the following line if loading pre-run benchmark data,
# will load into variable named "ks_bm"
# load(paste0(save_dir, "ks_bm.Rds"))

t_idx <- match("t_stop", colnames(ks_bm))
mbm_ks <- melt(ks_bm, measure.vars = -seq_len(t_idx),
               variable.name = "FuncName", value.name = "time")

Names <- c("ks_Gondan_14", "ks_Navarro_09")
Color <- c("#d7db42", "#33cc33")


### Violin Plot
ggplot(mbm_ks, aes(x = factor(FuncName, levels = Names), y = time,
                   color = factor(FuncName, levels = Names),
                   fill = factor(FuncName, levels = Names))) +
       geom_violin(trim = TRUE, alpha = 0.5) +
       scale_color_manual(values = Color, guide = FALSE) +
       scale_fill_manual(values = Color, guide = FALSE) +
       geom_boxplot(width = 0.15, fill = "white", alpha = 0.5) +
       stat_summary(fun = mean, geom = "errorbar",
                    aes(ymax = ..y.., ymin = ..y..),
                    width = .35, linetype = "dashed") +
       labs(title = "Distribution of median benchmark times",
            subtitle = "Dashed lines represent mean benchmark times",
            x = "Method", y = "Time (nanoseconds)",
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
ggsave("paper_analysis/images/ks_bm.png", width = 16, height = 9)


### Mean and Quantile Plots
# t
ggplot(mbm_ks, aes(x = t_mid, y = time,
                   color = factor(FuncName, levels = Names),
                   fill = factor(FuncName, levels = Names))) +
  stat_summary(fun.min = min, fun.max = max,
               geom = "ribbon", color = NA, alpha = 0.1) +
  stat_summary(fun.min = function(z) { quantile(z, 0.1) },
               fun.max = function(z) { quantile(z, 0.9) },
               geom = "ribbon", color = NA, alpha = 0.2) +
  stat_summary(fun = mean, geom = "line") +
  scale_color_manual(values = Color) +
  scale_fill_manual(values = Color) +
  labs(title = "Means of the median microbenchmark results",
       subtitle = paste(
         "The shaded regions represent the 10% and 90% quantiles",
         "The lighter shaded regions represent the min and max times",
         sep = ";\n"),
       x = bquote(frac(t, a^2) ~ ", effective response time"),
       y = "Time (microseconds)") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        plot.title = element_text(size = 20),
        plot.subtitle = element_text(size = 16,
                                     margin = margin(5, 5, 15, 5, "pt")),
        axis.text.x = element_text(size = 16, angle = 90,
                                   vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        strip.text = element_text(size = 14),
        strip.background = element_rect(fill = "white"),
        legend.position = "none") +
  facet_wrap(~ factor(FuncName, levels = Names))
ggsave("paper_analysis/images/ks_meq_t.png", width = 16, height = 9)
