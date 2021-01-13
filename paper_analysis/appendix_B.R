library("reshape2")
library("ggplot2")
library("ggforce")


########## violin plots, vector RT input
load("paper_analysis/benchmark_data/bm_vec.Rds")
t_idx <- match("SV", colnames(bm_vec))
bm_vec[, -seq_len(t_idx)] <- bm_vec[, -seq_len(t_idx)]/1000 # convert to microseconds
mbm_vec <- melt(bm_vec, measure.vars = -seq_len(t_idx),
                variable.name = "FuncName", value.name = "time")

Names_vec <- c("fb_SWSE_17", "fb_SWSE_14", "fb_Gon_17", "fb_Gon_14",
               "fb_Nav_17", "fb_Nav_14", "fs_SWSE_17", "fs_SWSE_14",
               "fs_Gon_17", "fs_Gon_14", "fs_Nav_17", "fs_Nav_14",
               "fl_Nav_09", "RWiener", "Gondan", "rtdists")
Color_vec <- c("#e000b4", "#ff99eb", "#e68a00", "#ffb366",
               "#006699", "#66ccff", "#9900cc", "#cc99ff",
               "#c2a500", "#d7db42", "#336600", "#33cc33",
               "#996633", "#ff9999", "#ff5050", "#990000")

mi <- min(bm_vec[, -seq_len(t_idx)])
ma <- max(bm_vec[, (t_idx+1):(ncol(bm_vec)-4)])

ggplot(mbm_vec, aes(x = factor(FuncName, levels = Names_vec), y = time,
                    color = factor(FuncName, levels = Names_vec),
                    fill = factor(FuncName, levels = Names_vec))) +
       geom_violin(trim = TRUE, alpha = 0.5) +
       scale_color_manual(values = Color_vec, guide = FALSE) +
       scale_fill_manual(values = Color_vec, guide = FALSE) +
       geom_boxplot(width = 0.15, fill = "white", alpha = 0.5) +
       stat_summary(fun = mean, geom = "errorbar",
                    aes(ymax = ..y.., ymin = ..y..),
                    width = .35, linetype = "dashed") +
       facet_zoom(ylim = c(mi, ma)) +
       labs(title = "Distribution of median benchmark times",
            subtitle = paste("Response time data input as a vector",
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
ggsave("paper_analysis/image_files/bm_vec_violin.png", width = 16, height = 9)



########## violin plots, individual RT input
load("paper_analysis/benchmark_data/bm_ind.Rds")
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
ggsave("paper_analysis/image_files/bm_ind_violin.png", width = 16, height = 9)



########## mean & quantile plots, individual RT input
Names_meq <- c("fb_SWSE_17", "fs_SWSE_14", "fl_Nav_09",
               "RWiener", "Gondan", "rtdists")
Color_meq <- c("#e000b4", "#cc99ff", "#996633",
               "#ff9999", "#ff5050", "#990000")
mbm_meq <- mbm_ind[mbm_ind[["FuncName"]] %in% Names_meq, ]


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
       subtitle = paste(
         "The shaded regions represent the 10% and 90% quantiles",
         "The lighter shaded regions represent the min and max times",
         sep = ";\n"),
       x = bquote(frac(rt, a^2) ~ ", effective response time, " ~ log[10]),
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
  facet_wrap(~ factor(FuncName, levels = Names_meq), scales = "free_y")
ggsave("paper_analysis/image_files/bm_ind_rtaa.png", width = 16, height = 9)


ggplot(mbm_meq, aes(x = W, y = time,
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
  labs(title = "Means of the median microbenchmark results",
       subtitle = paste(
         "The shaded regions represent the 10% and 90% quantiles",
         "The lighter shaded regions represent the min and max times",
         sep = ";\n"),
       x = "w, relative starting point (a priori bias)",
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
  facet_wrap(~ factor(FuncName, levels = Names_meq), scales = "free_y")
ggsave("paper_analysis/image_files/bm_ind_w.png", width = 16, height = 9)


ggplot(mbm_meq, aes(x = V, y = time,
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
  labs(title = "Means of the median microbenchmark results",
       subtitle = paste(
         "The shaded regions represent the 10% and 90% quantiles",
         "The lighter shaded regions represent the min and max times",
         sep = ";\n"),
       x = "v, drift rate",
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
  facet_wrap(~ factor(FuncName, levels = Names_meq), scales = "free_y")
ggsave("paper_analysis/image_files/bm_ind_v.png", width = 16, height = 9)


ggplot(mbm_meq, aes(x = SV, y = time,
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
  labs(title = "Means of the median microbenchmark results",
       subtitle = paste(
         "The shaded regions represent the 10% and 90% quantiles",
         "The lighter shaded regions represent the min and max times",
         sep = ";\n"),
       x = "sv, inter-trial variability in the drift rate",
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
  facet_wrap(~ factor(FuncName, levels = Names_meq), scales = "free_y")
ggsave("paper_analysis/image_files/bm_ind_sv.png", width = 16, height = 9)






########## violin plots, individual RT input (0-2 seconds from Section 4.2)
load("paper_analysis/benchmark_data/bm_ind_0-2.Rds")
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
ggsave("paper_analysis/image_files/bm_ind_0-2_violin.png", width = 16, height = 9)
