library("ggplot2")


# Import fitting data
fit <- readRDS("benchmark_testing/Results/fit.Rds")
fit$BmTime <- fit$BmTime * 1e-6 # convert to seconds

# Prepare vectors of algorithm names and associated colors
Names <- c("fs_Fos_17", "fs_Kes_17", "fs_Nav_17",
           "fb_Kes_17", "fb_Nav_17",
           "RWiener", "Kesselmeier", "rtdists")
Color <- c("#9900cc", "#006699", "#336600",
           "#c2a500", "#e68a00",
           "#ff9999", "#ff5050", "#990000")


# Plot number of calls for each algorithm
calls_by_alg <- ggplot(fit, aes(x = Algorithm, y = AlgoCalls,
                                color = factor(Algorithm, levels = Names),
                                fill = factor(Algorithm, levels = Names))) +
                geom_violin(trim = TRUE) +
                scale_color_manual(values = Color) +
                scale_fill_manual(values = Color) +
                geom_boxplot(width = 0.15, fill="white") +
                labs(title = "Number of calls to algorithms during fitting",
                     x = "Algorithm", y = "Algorithm Calls",
                     color = "Algorithm") +
                theme_bw() +
                theme(panel.grid.minor = element_blank(),
                      panel.border = element_blank(),
                      plot.title = element_text(size = 23),
                      axis.text.x = element_text(size = 16, angle = 30),
                      axis.text.y = element_text(size = 16),
                      axis.title.x = element_text(size = 20),
                      axis.title.y = element_text(size = 20),
                      legend.position = "none")
ggsave("benchmark_testing/Results/Images/calls_by_alg.png", plot = calls_by_alg,
       width = 16, height = 9)


#  Mean benchmark times (violin plots) for each method
violin <- ggplot(fit, aes(x = Algorithm, y = BmTime,
                          color = factor(Algorithm, levels = Names),
                          fill = factor(Algorithm, levels = Names))) +
                 geom_violin(trim = TRUE, alpha = 0.5) +
                 scale_color_manual(values = Color) +
                 scale_fill_manual(values = Color) +
                 geom_boxplot(width = 0.15, fill = "white", alpha = 0.5) +
                 stat_summary(fun.y = mean, geom = "errorbar",
                              aes(ymax = ..y.., ymin = ..y..),
                              width = .35, linetype = "dashed") +
                 # abline
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
ggsave("benchmark_testing/Results/Images/bm_by_alg.png", plot = violin,
       width = 16, height = 9)
