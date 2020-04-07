library("reshape2")
library("ggplot2")
library("tidyverse")




############## Simulated Data ##############
### Complementary ECDF Plots
bm <- readRDS("benchmark_testing/Results/sim_vec_10000.Rds")
t_idx <- match("FuncName", colnames(bm))
temp <- bm[,-seq_len(t_idx)]
bm[,-seq_len(t_idx)] <- bm[, -seq_len(t_idx)]/1000
mbm <- melt(bm, measure.vars = -seq_len(t_idx),
            variable.name="iter", value.name="time")


Names_RCp <- c("fs_Fos_17", "fs_Fos_14", "fs_Kes_17", "fs_Kes_14",
               "fs_Nav_17", "fs_Nav_14", "fl_Nav_09",
               "fb_Kes_17", "fb_Kes_14", "fb_Nav_17", "fb_Nav_14",
               "rtdists", "RWiener", "Kesselmeier")
Color_RCp <- c("#9900cc", "#cc99ff", "#006699", "#66ccff",
               "#336600", "#33cc33", "#996633",
               "#c2a500", "#d7db42", "#e68a00", "#ffb366",
               "#990000", "#ff9999", "#ff5050")
mbm_RCp <- mbm
vec_ecdf_RCp <- ggplot(mbm_RCp, aes(x=time, color=factor(FuncName, levels=Names_RCp))) +
  geom_line(aes(y = 1 - ..y..), stat='ecdf') +
  scale_color_manual(values=Color_RCp) +
  coord_cartesian(xlim=c(25,400)) +
  labs(title="Complementary ECDF of Microbenchmark Results",
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
        legend.text = element_text(size = 16, angle=0))
ggsave("benchmark_testing/Results/Images/vec_10000_ecdf_RCp.png", plot=vec_ecdf_RCp, width=16, height=9)


Names_sum <- c("fs_Fos_17", "fs_Fos_14", "fs_Kes_17", "fs_Kes_14",
               "fs_Nav_17", "fs_Nav_14",
               "fb_Kes_17", "fb_Kes_14", "fb_Nav_17", "fb_Nav_14")
Color_sum <- c("#9900cc", "#cc99ff", "#006699", "#66ccff",
               "#336600", "#33cc33",
               "#c2a500", "#d7db42", "#e68a00", "#ffb366")
mbm_sum <- subset(mbm, FuncName %in% Names_sum)
vec_ecdf_sum <- ggplot(mbm_sum, aes(x=time, color=factor(FuncName, levels=Names_sum))) +
  geom_line(aes(y = 1 - ..y..), stat='ecdf') +
  scale_color_manual(values=Color_sum) +
  coord_cartesian(xlim=c(25,200)) +
  labs(title="Complementary ECDF of Microbenchmark Results",
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
        legend.text = element_text(size = 16, angle=0))
ggsave("benchmark_testing/Results/Images/vec_10000_ecdf_sum.png", plot=vec_ecdf_sum, width=16, height=9)


Names_slt <- c("fs_Fos_17", "fs_Kes_17",
               "fs_Nav_17", "fl_Nav_09",
               "fb_Kes_17", "fb_Nav_17")
Color_slt <- c("#9900cc", "#006699",
               "#336600", "#996633",
               "#c2a500", "#e68a00")
mbm_slt <- subset(mbm, FuncName %in% Names_slt)
vec_ecdf_slt <- ggplot(mbm_slt, aes(x=time, color=factor(FuncName, levels=Names_slt))) +
  geom_line(aes(y = 1 - ..y..), stat='ecdf') +
  scale_color_manual(values=Color_slt) +
  coord_cartesian(xlim=c(25,200)) +
  labs(title="Complementary ECDF of Microbenchmark Results",
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
        legend.text = element_text(size = 16, angle=0))
ggsave("benchmark_testing/Results/Images/vec_10000_ecdf_slt.png", plot=vec_ecdf_slt, width=16, height=9)


Names_bRC <- c("fs_Fos_17", "fb_Kes_17", "fb_Nav_17",
               "rtdists", "RWiener", "Kesselmeier")
Color_bRC <- c("#9900cc", "#c2a500", "#e68a00",
               "#990000", "#ff9999", "#ff5050")
mbm_bRC <- subset(mbm, FuncName %in% Names_bRC)
vec_ecdf_bRC <- ggplot(mbm_bRC, aes(x=time, color=factor(FuncName, levels=Names_bRC))) +
  geom_line(aes(y = 1 - ..y..), stat='ecdf') +
  scale_color_manual(values=Color_bRC) +
  coord_cartesian(xlim=c(25,200)) +
  labs(title="Complementary ECDF of Microbenchmark Results",
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
        legend.text = element_text(size = 16, angle=0))
ggsave("benchmark_testing/Results/Images/vec_10000_ecdf_bRC.png", plot=vec_ecdf_bRC, width=16, height=9)













### Each parameter vs benchmark time
bm <- readRDS("benchmark_testing/Results/sim_10000.Rds")
t_idx <- match("FuncName", colnames(bm))
temp <- bm[,-seq_len(t_idx)]
bm[,-seq_len(t_idx)] <- bm[, -seq_len(t_idx)]/1000
mbm <- melt(bm, measure.vars = -seq_len(t_idx),
            variable.name="iter", value.name="time")

Names_bRC <- c("fs_Fos_17", "fb_Kes_17", "fb_Nav_17",
               "rtdists", "RWiener", "Kesselmeier")
Color_bRC <- c("#9900cc", "#c2a500", "#e68a00",
               "#990000", "#ff9999", "#ff5050")
mbm_bRC <- subset(mbm, FuncName %in% Names_bRC)

# RT
mbm_bRC %>%
  mutate(FuncName = factor(FuncName, levels=Names_bRC)) %>%
  group_by(FuncName, RT) %>%
  summarise(means = mean(time),
            upper = quantile(time, prob = 0.9),
            lower = quantile(time, prob = 0.1)) %>%
  ggplot(aes(x = RT, y = means, color = FuncName)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = FuncName),
              alpha = 0.2, color = NA) +
  geom_line(aes(group = 1)) +
  scale_color_manual(values = Color_bRC) +
  scale_fill_manual(values = Color_bRC) +
  labs(title = "Means with 10% and 90% Quantiles of Microbenchmark Results",
       x = "RT, Response Time", y = "Time (ms)") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        strip.background = element_rect(fill = "white"),
        legend.position = "none") +
  facet_wrap("FuncName")
ggsave("benchmark_testing/Results/Images/meq_RT.png", width=16, height=9)

# A
mbm_bRC %>%
  mutate(FuncName = factor(FuncName, levels=Names_bRC)) %>%
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
  scale_color_manual(values = Color_bRC) +
  scale_fill_manual(values = Color_bRC) +
  labs(title = "Means with 10% and 90% Quantiles of Microbenchmark Results",
       subtitle = "Grey shaded region is typical range of the parameter",
       x = "A, Threshold Separation", y = "Time (ms)") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        strip.background = element_rect(fill = "white"),
        legend.position = "none") +
  facet_wrap("FuncName")
ggsave("benchmark_testing/Results/Images/meq_A.png", width=16, height=9)

# V
mbm_bRC %>%
  mutate(FuncName = factor(FuncName, levels=Names_bRC)) %>%
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
  scale_color_manual(values = Color_bRC) +
  scale_fill_manual(values = Color_bRC) +
  labs(title = "Means with 10% and 90% Quantiles of Microbenchmark Results",
       subtitle = "Grey shaded region is typical range of the parameter",
       x = "V, Drift Rate", y = "Time (ms)") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        strip.background = element_rect(fill = "white"),
        legend.position = "none") +
  facet_wrap("FuncName")
ggsave("benchmark_testing/Results/Images/meq_V.png", width=16, height=9)

# W
mbm_bRC %>%
  mutate(FuncName = factor(FuncName, levels=Names_bRC)) %>%
  group_by(FuncName, W) %>%
  summarise(means = mean(time),
            upper = quantile(time, prob = 0.9),
            lower = quantile(time, prob = 0.1)) %>%
  ggplot(aes(x = W, y = means, color = FuncName)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = FuncName),
              alpha = 0.2, color = NA) +
  geom_line(aes(group = 1)) +
  scale_color_manual(values = Color_bRC) +
  scale_fill_manual(values = Color_bRC) +
  labs(title = "Means with 10% and 90% Quantiles of Microbenchmark Results",
       x = "W, Relative Starting Point (A Priori Bias)", y = "Time (ms)") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        strip.background = element_rect(fill = "white"),
        legend.position = "none") +
  facet_wrap("FuncName")
ggsave("benchmark_testing/Results/Images/meq_W.png", width=16, height=9)













### Mean benchmark times (violin plots) for each method
bm <- readRDS("benchmark_testing/Results/sim_vec_10000.Rds")
t_idx <- match("FuncName", colnames(bm))
temp <- bm[,-seq_len(t_idx)]
bm[,-seq_len(t_idx)] <- bm[, -seq_len(t_idx)]/1000
mbm <- melt(bm, measure.vars = -seq_len(t_idx),
            variable.name="iter", value.name="time")


Names <- c("fs_Fos_17", "fs_Fos_14", "fs_Kes_17", "fs_Kes_14",
           "fs_Nav_17", "fs_Nav_14", "fl_Nav_09",
           "fb_Kes_17", "fb_Kes_14", "fb_Nav_17", "fb_Nav_14",
           "rtdists", "RWiener", "Kesselmeier")
Color <- c("#9900cc", "#cc99ff", "#006699", "#66ccff",
           "#336600", "#33cc33", "#996633",
           "#c2a500", "#d7db42", "#e68a00", "#ffb366",
           "#990000", "#ff9999", "#ff5050")

sim_vec_viol <- ggplot(mbm, aes(x = FuncName, y = time,
                                color = factor(FuncName, levels = Names),
                                fill = factor(FuncName, levels = Names))) +
                geom_violin(trim = TRUE) +
                scale_color_manual(values = Color) +
                scale_fill_manual(values = Color) +
                geom_boxplot(width = 0.15, fill="white") +
                coord_cartesian(ylim=c(0,600)) +
                labs(title = "Distribution of benchmark times",
                     x = "Method", y = "Time (ms)",
                     color = "Method") +
                theme_bw() +
                theme(panel.grid.minor = element_blank(),
                      panel.border = element_blank(),
                      plot.title = element_text(size = 23),
                      axis.text.x = element_text(size = 16, angle = 30),
                      axis.text.y = element_text(size = 16),
                      axis.title.x = element_text(size = 20),
                      axis.title.y = element_text(size = 20),
                      legend.position = "none")
ggsave("benchmark_testing/Results/Images/sim_vec_viol.png", plot=sim_vec_viol,
       width=16, height=9)
