library("ggplot2")


# Import fitting data
fit <- readRDS("benchmark_testing/Results/fit.Rds")


# Prepare vectors of algorithm names and associated colors
Names <- c("fs_Fos_17", "fb_Kes_17", "fb_Nav_17",
           "rtdists", "RWiener", "Kesselmeier")
Color <- c("#9900cc", "#c2a500", "#e68a00",
           "#990000", "#ff9999", "#ff5050")


# Plot number of calls for each algorithm
calls_by_alg <- ggplot(fit, aes(x = Algorithm, y = AlgoCalls,
                                color = factor(Algorithm, levels = Names),
                                fill = factor(Algorithm, levels = Names))) +
                geom_violin(trim = FALSE) +
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
                      axis.text.x = element_text(size = 16),
                      axis.text.y = element_text(size = 16),
                      axis.title.x = element_text(size = 20),
                      axis.title.y = element_text(size = 20),
                      legend.position = "none")
ggsave("benchmark_testing/Results/Images/calls_by_alg.png", plot=calls_by_alg,
       width=16, height=9)


# Plot number of calls (for each algorithm) for each individual
calls_by_ind <- ggplot(fit, aes(x = ind, y = AlgoCalls,
                                color = factor(Algorithm, levels = Names))) +
                geom_point(alpha = 0.2) +
                scale_color_manual(values=Color) +
                labs(title = "Number of calls to algorithms during fitting, for each individual",
                     x = "Individual", y = "Algorithm Calls",
                     color = "Algorithm") +
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
                      legend.text = element_text(size = 16, angle = 0))
ggsave("benchmark_testing/Results/Images/calls_by_ind.png", plot = calls_by_ind,
       width = 16, height = 9)
