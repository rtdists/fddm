library("reshape2")
library("ggplot2")
library("tidyverse")




############## Simulated Data ##############

### Mean benchmark times (violin plots) for each method
bm <- readRDS("benchmark_testing/Results/sim_0-10_10000.Rds")
t_idx <- match("W", colnames(bm))
bm[,-seq_len(t_idx)] <- bm[, -seq_len(t_idx)]/1000
mbm <- melt(bm, measure.vars = -seq_len(t_idx),
            variable.name = "FuncName", value.name = "time")


Names <- c("fddm_fast", "fs_Fos_17", "fs_Fos_14",
           "fs_Kes_17", "fs_Kes_14", "fs_Nav_17", "fs_Nav_14",
           "fb_Kes_17", "fb_Kes_14", "fb_Nav_17", "fb_Nav_14",
           "fl_Nav_09", "RWiener", "Kesselmeier", "rtdists")
Color <- c("#ff00cc", "#9900cc", "#cc99ff",
           "#006699", "#66ccff", "#336600", "#33cc33",
           "#c2a500", "#d7db42", "#e68a00", "#ffb366",
           "#996633", "#ff9999", "#ff5050", "#990000")

violin <- ggplot(mbm, aes(x = FuncName, y = time,
                              color = factor(FuncName, levels = Names),
                              fill = factor(FuncName, levels = Names))) +
                  geom_violin(trim = TRUE, alpha = 0.5) +
                  scale_color_manual(values = Color) +
                  scale_fill_manual(values = Color) +
                  geom_boxplot(width = 0.15, fill = "white", alpha = 0.5) +
                  stat_summary(fun.y = mean, geom = "errorbar",
                               aes(ymax = ..y.., ymin = ..y..),
                               width = .35, linetype = "dashed") +
                  # abline
                  coord_cartesian(ylim = c(0,50)) +
                  labs(title = "Distribution of median benchmark times",
                       subtitle = "Dashed lines represent mean benchmark times",
                       x = "Method", y = "Time (ms)",
                       color = "Method") +
                  theme_bw() +
                  theme(#panel.grid.minor = element_blank(),
                        panel.border = element_blank(),
                        plot.title = element_text(size = 23),
                        plot.subtitle = element_text(size = 16),
                        axis.text.x = element_text(size = 16, angle = 30),
                        axis.text.y = element_text(size = 16),
                        axis.title.x = element_text(size = 20),
                        axis.title.y = element_text(size = 20),
                        legend.position = "none")
ggsave("benchmark_testing/Results/Images/violin_0-10_indRT.png", plot = violin,
       width = 16, height = 9)



### Densities of parameters in the tails of the violins (for each method)
plot_par_dens <- function(df, fname) {
  # choose the method
  fbm <- df[, c("V", "A", "W", fname)]
  times <- fbm[,ncol(fbm)]
  # crudely extract tail
  tail_start <- median(times) + 2*sd(times)
  ftbm <- (fbm[fbm[,ncol(fbm)] > tail_start,])[,-ncol(fbm)]
  # plot density of parameters
  ftbm_melt <- melt(ftbm, id.vars = NULL, variable.name = "Parameter")
  ftbm_par_dens <- ggplot(ftbm_melt, aes(value)) +
                          geom_density() +
                          labs(title = "Tail densities of parameters",
                               subtitle = fname,
                               x = "Parameter Value", y = "Density") +
                          theme_bw() +
                          theme(panel.grid.minor = element_blank(),
                                panel.border = element_blank(),
                                strip.background = element_rect(fill = "white"),
                                legend.position = "none") +
                          facet_wrap("Parameter", scales = "free")
  return(ftbm_par_dens)
}
plot_par_dens(bm, Names[2]) # Names[2:15]
ggsave("benchmark_testing/Results/Images/violin_tails/fs_Fos_17_0-10.png",
       width = 16, height = 9)









### Each parameter vs benchmark time
bm <- readRDS("benchmark_testing/Results/sim_0-10_10000.Rds")
t_idx <- match("W", colnames(bm))
bm[,-seq_len(t_idx)] <- bm[, -seq_len(t_idx)]/1000
mbm <- melt(bm, measure.vars = -seq_len(t_idx),
            variable.name = "FuncName", value.name = "time")

Names_meq <- c("fs_Fos_17", "fs_Kes_17", "fs_Nav_17",
               "fddm_fast", "fb_Kes_17", "fb_Nav_17",
               "RWiener", "Kesselmeier", "fl_Nav_09")
Color_meq <- c("#9900cc", "#006699", "#336600",
               "#ff00cc", "#c2a500", "#e68a00",
               "#ff9999", "#ff5050", "#996633")
mbm_meq <- subset(mbm, FuncName %in% Names_meq)

# RT
mbm_meq %>%
  mutate(FuncName = factor(FuncName, levels = Names_meq)) %>%
  group_by(FuncName, RT) %>%
  summarise(means = mean(time),
            upper = quantile(time, prob = 0.9),
            lower = quantile(time, prob = 0.1),
            max = max(time),
            min = min(time)) %>%
  ggplot(aes(x = RT, y = means, color = FuncName)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = FuncName),
              alpha = 0.2, color = NA) +
  geom_ribbon(aes(ymin = min, ymax = max, fill = FuncName),
              alpha = 0.1, color = NA) +
  geom_line(aes(group = 1)) +
  scale_color_manual(values = Color_meq) +
  scale_fill_manual(values = Color_meq) +
  labs(title = "Means with 10% and 90% quantiles of median microbenchmark results",
       x = "rt, response time", y = "Time (ms)") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        strip.background = element_rect(fill = "white"),
        legend.position = "none") +
  facet_wrap("FuncName")
ggsave("benchmark_testing/Results/Images/meq_plots/meq_RT_0-10.png",
       width = 16, height = 9)

# A
mbm_meq %>%
  mutate(FuncName = factor(FuncName, levels = Names_meq)) %>%
  group_by(FuncName, A) %>%
  summarise(means = mean(time),
            upper = quantile(time, prob = 0.9),
            lower = quantile(time, prob = 0.1)) %>%
  ggplot(aes(x = A, y = means, color = FuncName)) +
  annotate("rect", fill = "lightgrey", alpha = 0.3,
           xmin = 0.5, xmax = 2, ymin = -Inf, ymax = Inf) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = FuncName),
              alpha = 0.2, color = NA) +
  geom_line(aes(group = 1)) +
  scale_color_manual(values = Color_meq) +
  scale_fill_manual(values = Color_meq) +
  labs(title = "Means with 10% and 90% quantiles of median microbenchmark results",
       subtitle = "Grey shaded region is typical range of the parameter",
       x = "a, threshold separation", y = "Time (ms)") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        strip.background = element_rect(fill = "white"),
        legend.position = "none") +
  facet_wrap("FuncName")
ggsave("benchmark_testing/Results/Images/meq_plots/meq_A_0-10.png",
       width = 16, height = 9)

# V
mbm_meq %>%
  mutate(FuncName = factor(FuncName, levels = Names_meq)) %>%
  group_by(FuncName, V) %>%
  summarise(means = mean(time),
            upper = quantile(time, prob = 0.9),
            lower = quantile(time, prob = 0.1)) %>%
  ggplot(aes(x = V, y = means, color = FuncName)) +
  annotate("rect", fill = "lightgrey", alpha = 0.3,
           xmin = -5, xmax = 5, ymin = -Inf, ymax = Inf) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = FuncName),
              alpha = 0.2, color = NA) +
  geom_line(aes(group = 1)) +
  scale_color_manual(values = Color_meq) +
  scale_fill_manual(values = Color_meq) +
  labs(title = "Means with 10% and 90% quantiles of median microbenchmark results",
       subtitle = "Grey shaded region is typical range of the parameter",
       x = "v, drift rate", y = "Time (ms)") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        strip.background = element_rect(fill = "white"),
        legend.position = "none") +
  facet_wrap("FuncName")
ggsave("benchmark_testing/Results/Images/meq_plots/meq_V_0-10.png",
       width = 16, height = 9)

# W
mbm_meq %>%
  mutate(FuncName = factor(FuncName, levels = Names_meq)) %>%
  group_by(FuncName, W) %>%
  summarise(means = mean(time),
            upper = quantile(time, prob = 0.9),
            lower = quantile(time, prob = 0.1)) %>%
  ggplot(aes(x = W, y = means, color = FuncName)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = FuncName),
              alpha = 0.2, color = NA) +
  geom_line(aes(group = 1)) +
  scale_color_manual(values = Color_meq) +
  scale_fill_manual(values = Color_meq) +
  labs(title = "Means with 10% and 90% quantiles of median microbenchmark results",
       x = "w, relative starting point (a priori bias)", y = "Time (ms)") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        strip.background = element_rect(fill = "white"),
        legend.position = "none") +
  facet_wrap("FuncName")
ggsave("benchmark_testing/Results/Images/meq_plots/meq_W_0-10.png",
       width = 16, height = 9)










### Complementary ECDF Plots
bm <- readRDS("benchmark_testing/Results/vec_0-3_1000.Rds")
t_idx <- match("W", colnames(bm))
bm[,-seq_len(t_idx)] <- bm[, -seq_len(t_idx)]/1000
mbm <- melt(bm, measure.vars = -seq_len(t_idx),
            variable.name = "FuncName", value.name = "time")


Names_RCp <- c("fddm_fast", "fs_Fos_17", "fs_Fos_14",
               "fs_Kes_17", "fs_Kes_14", "fs_Nav_17", "fs_Nav_14",
               "fb_Kes_17", "fb_Kes_14", "fb_Nav_17", "fb_Nav_14",
               "fl_Nav_09", "RWiener", "Kesselmeier", "rtdists")
Color_RCp <- c("#ff00cc", "#9900cc", "#cc99ff",
               "#006699", "#66ccff", "#336600", "#33cc33",
               "#c2a500", "#d7db42", "#e68a00", "#ffb366",
               "#996633", "#ff9999", "#ff5050", "#990000")
mbm_RCp <- mbm
ecdf_RCp <- ggplot(mbm_RCp, aes(x = time,
                                color = factor(FuncName, levels = Names_RCp))) +
  geom_line(aes(y = 1 - ..y..), stat = 'ecdf') +
  scale_color_manual(values=Color_RCp) +
  # coord_cartesian(xlim = c(0,50)) +
  labs(title = "Complementary ECDF of median Microbenchmark results",
       x = "Time (ms)", y = "Complementary ECDF",
       color = "Method") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        plot.title = element_text(size = 23),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.position = c(1,1),
        legend.justification = c(1,1),
        legend.direction = "vertical",
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 16, angle = 0)) +
  guides(color = guide_legend(override.aes = list(size = 2)))
ggsave("benchmark_testing/Results/Images/ecdf_plots/ecdf_RCp_0-3.png",
       plot = ecdf_RCp, width = 16, height = 9)


Names_sum <- c("fs_Fos_17", "fs_Fos_14",
               "fs_Kes_17", "fs_Kes_14", "fs_Nav_17", "fs_Nav_14",
               "fb_Kes_17", "fb_Kes_14", "fb_Nav_17", "fb_Nav_14")
Color_sum <- c("#9900cc", "#cc99ff",
               "#006699", "#66ccff", "#336600", "#33cc33",
               "#c2a500", "#d7db42", "#e68a00", "#ffb366")
mbm_sum <- subset(mbm, FuncName %in% Names_sum)
ecdf_sum <- ggplot(mbm_sum, aes(x = time,
                                color = factor(FuncName, levels = Names_sum))) +
  geom_line(aes(y = 1 - ..y..), stat = 'ecdf') +
  scale_color_manual(values = Color_sum) +
  # coord_cartesian(xlim = c(25,200)) +
  labs(title = "Complementary ECDF of median Microbenchmark results",
       x = "Time (ms)", y = "Complementary ECDF",
       color = "Method") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        plot.title = element_text(size = 23),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.position = c(1,1),
        legend.justification = c(1,1),
        legend.direction = "vertical",
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 16, angle = 0)) +
  guides(color = guide_legend(override.aes = list(size = 2)))
ggsave("benchmark_testing/Results/Images/ecdf_plots/ecdf_sum_0-3.png",
       plot = ecdf_sum, width = 16, height = 9)


Names_slt <- c("fs_Fos_17", "fs_Kes_17", "fs_Nav_17",
               "fb_Kes_17", "fb_Nav_17", "fl_Nav_09")
Color_slt <- c("#9900cc", "#006699", "#336600",
               "#c2a500", "#e68a00", "#996633")
mbm_slt <- subset(mbm, FuncName %in% Names_slt)
ecdf_slt <- ggplot(mbm_slt, aes(x = time,
                                color = factor(FuncName, levels = Names_slt))) +
  geom_line(aes(y = 1 - ..y..), stat = 'ecdf') +
  scale_color_manual(values = Color_slt) +
  # coord_cartesian(xlim = c(25,200)) +
  labs(title = "Complementary ECDF of median Microbenchmark results",
       x = "Time (ms)", y = "Complementary ECDF",
       color="Method") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        plot.title = element_text(size = 23),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.position = c(1,1),
        legend.justification = c(1,1),
        legend.direction = "vertical",
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 16, angle = 0)) +
  guides(color = guide_legend(override.aes = list(size = 2)))
ggsave("benchmark_testing/Results/Images/ecdf_plots/ecdf_slt_0-3.png",
       plot = ecdf_slt, width = 16, height = 9)


Names_bRC <- c("fddm_fast", "fs_Fos_17", "fs_Kes_17", "fs_Nav_17",
               "fb_Kes_17", "fb_Nav_17", "RWiener", "Kesselmeier")
Color_bRC <- c("#ff00cc", "#9900cc", "#006699", "#336600",
               "#c2a500", "#e68a00", "#ff9999", "#ff5050")
mbm_bRC <- subset(mbm, FuncName %in% Names_bRC)
ecdf_bRC <- ggplot(mbm_bRC, aes(x = time,
                                color = factor(FuncName, levels = Names_bRC))) +
  geom_line(aes(y = 1 - ..y..), stat = 'ecdf') +
  scale_color_manual(values = Color_bRC) +
  # coord_cartesian(xlim = c(25,200)) +
  labs(title="Complementary ECDF of median Microbenchmark results",
       x="Time (ms)", y="Complementary ECDF",
       color="Method") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        plot.title = element_text(size = 23),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.position = c(1,1),
        legend.justification = c(1,1),
        legend.direction = "vertical",
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 16, angle = 0)) +
  guides(color = guide_legend(override.aes = list(size = 2)))
ggsave("benchmark_testing/Results/Images/ecdf_plots/ecdf_bRC_0-3_TEST2.png",
       plot = ecdf_bRC, width = 16, height = 9)
