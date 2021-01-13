library("reshape2")
library("ggplot2")
library("ggforce")


########## violin plots, vector RT input (0-2 seconds)
load("paper_analysis/benchmark_data/bm_vec_0-2.Rds")
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
ggsave("paper_analysis/image_files/bm_vec_0-2_violin.png", width = 16, height = 9)



########## violin plots, individual RT input (0-2 seconds also in appendix_B.R)
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
ggsave("paper_analysis/image_files/bm_ind_rtaa_0-2.png", width = 16, height = 9)


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
ggsave("paper_analysis/image_files/bm_ind_w_0-2.png", width = 16, height = 9)


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
ggsave("paper_analysis/image_files/bm_ind_v_0-2.png", width = 16, height = 9)


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
ggsave("paper_analysis/image_files/bm_ind_sv_0-2.png", width = 16, height = 9)






########## fitting benchmark plots
fit_prep <- function(fit, eps = 1e-4) {
  nr <- nrow(fit)
  fit[["Obj_diff"]] <- rep(0, nr)

  ids <- unique(fit[["ID"]])
  nids <- length(ids)
  algos <- unique(fit[["Algorithm"]])
  nalgos <- length(algos)

  ninit <- nrow(fit[fit[["ID"]] == ids[1] & fit[["Algorithm"]] == algos[1], ])
  for (i in 1:nids) {
    for (j in 1:ninit) {
      idx <- which(fit[["ID"]] == ids[i])[ninit*(0:(nalgos-1)) + j]
      objs <- fit[idx, "Objective"]
      min_obj <- min(objs)
      abs_min_obj <- abs(min_obj)
      obj_diffs <- objs - min(objs)
      fit[idx, "Obj_diff"] <- ifelse(obj_diffs <= eps*abs_min_obj, 0,
        ifelse(obj_diffs > eps*abs_min_obj & obj_diffs <= 2*abs_min_obj, 1, 3))
    }
  }

  fit[["BmTime"]] <- fit[["BmTime"]]*1e-6 # convert to milliseconds
  fit[["Convergence"]] <- ifelse(fit[["Convergence"]] < 1, 0, 1)

  return(fit)
}

# load data, will be in the variable 'fit'
load("paper_analysis/benchmark_data/bm_fit.Rds")
fit <- fit_prep(fit)

Names <- c("fb_SWSE_17", "fb_Gon_17", "fs_SWSE_14",
           "fs_Gon_17", "fl_Nav_09", "rtdists")
Color <- c("#e000b4", "#e68a00", "#cc99ff",
          "#c2a500", "#996633", "#990000")
Shape <- c(21, 25)
Sizes <- c(0, 3, 3)
Stroke <- c(0, 1, 1)
Fills <- c("#ffffff00", "#ffffff00", "#80808099")


fit_mbm <- melt(fit, id.vars = c("Algorithm", "Convergence", "Obj_diff"),
                measure.vars = "BmTime", value.name = "BmTime")[,-4]

mibm <- min(fit[fit[["Algorithm"]] != "rtdists", "BmTime"])
mabm <- max(fit[fit[["Algorithm"]] != "rtdists", "BmTime"])

ggplot(fit_mbm, aes(x = factor(Algorithm, levels = Names),
                    y = BmTime)) +
  geom_violin(trim = TRUE, alpha = 0.5,
              aes(color = factor(Algorithm, levels = Names),
                  fill = factor(Algorithm, levels = Names))) +
  geom_boxplot(width = 0.2, outlier.shape = NA,
               fill = "white", alpha = 0.4,
               aes(color = factor(Algorithm, levels = Names))) +
  stat_summary(fun = mean, geom = "errorbar",
               aes(ymax = ..y.., ymin = ..y..),
               width = .5, linetype = "dashed",
               color = Color) +
  scale_color_manual(values = Color, guide = FALSE) +
  scale_shape_manual(values = Shape,
                     name = "Convergence Code",
                     breaks = c(0, 1),
                     labels = c("Success", "Failure")) +
  scale_size_manual(values = Sizes, guide = FALSE) +
  scale_discrete_manual(aesthetics = "stroke", values = Stroke, guide = FALSE) +
  scale_fill_manual(values = Color, guide = FALSE) +
  ggnewscale::new_scale_fill() +
  scale_fill_manual(values = Fills,
                    name = "Objective Difference",
                    breaks = c(1, 2, 3),
                    labels = c("< 2", "NA", "> 2")) +
  geom_point(aes(color = factor(Algorithm, levels = Names),
                 shape = factor(Convergence, levels = c(0, 1)),
                 size = factor(Obj_diff, levels = c(0, 1, 3)),
                 stroke = factor(Obj_diff, levels = c(0, 1, 3)),
                 fill = factor(Obj_diff, levels = c(0, 1, 3)))) +
  coord_cartesian(ylim = c(mibm, mabm)) +
  labs(title = "Median microbenchmark times for data fitting",
       subtitle = paste(
         "Dashed lines represent mean benchmark times",
         "Visible points have an objective difference greater than 1e-4",
         sep = ";\n"),
       x = "Method", y = "Time (milliseconds)") +
  guides(shape = guide_legend(order = 1,
                              override.aes = list(size = Sizes[c(2, 3)])),
         fill = guide_legend(order = 2,
                             override.aes = list(size = Sizes[c(2, 3)],
                                                 shape = c(21, 21),
                                                 fill = Fills[c(2, 3)]))) +
  theme_bw() +
  theme(panel.border = element_blank(),
        plot.title = element_text(size = 23),
        plot.subtitle = element_text(size = 16,
                                     margin = margin(5, 5, 10, 5, "pt")),
        axis.text.x = element_text(size = 16, angle = 90,
                                   vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20,
                                    margin = margin(5, 10, 5, 5, "pt")),
        legend.position = "right",
        legend.box = "vertical",
        legend.direction = "vertical",
        legend.background = element_rect(fill = "transparent"),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 13))
ggsave("paper_analysis/image_files/bm_fit_bm.png", width = 16, height = 9)



fit_fev <- melt(fit, id.vars = c("Algorithm", "Convergence", "Obj_diff"),
                measure.vars = "FuncEvals", value.name = "FuncEvals")[,-4]

mife <- min(fit[fit[["Algorithm"]] != "rtdists", "FuncEvals"])
mafe <- max(fit[fit[["Algorithm"]] != "rtdists", "FuncEvals"])

ggplot(fit_fev, aes(x = factor(Algorithm, levels = Names),
                    y = FuncEvals)) +
  geom_violin(trim = TRUE, alpha = 0.5,
              aes(color = factor(Algorithm, levels = Names),
                  fill = factor(Algorithm, levels = Names))) +
  geom_boxplot(width = 0.2, outlier.shape = NA,
               fill = "white", alpha = 0.4,
               aes(color = factor(Algorithm, levels = Names))) +
  stat_summary(fun = mean, geom = "errorbar",
               aes(ymax = ..y.., ymin = ..y..),
               width = .5, linetype = "dashed",
               color = Color) +
  scale_color_manual(values = Color, guide = FALSE) +
  scale_shape_manual(values = Shape,
                     name = "Convergence Code",
                     breaks = c(0, 1),
                     labels = c("Success", "Failure")) +
  scale_size_manual(values = Sizes, guide = FALSE) +
  scale_discrete_manual(aesthetics = "stroke", values = Stroke, guide = FALSE) +
  scale_fill_manual(values = Color, guide = FALSE) +
  ggnewscale::new_scale_fill() +
  scale_fill_manual(values = Fills,
                    name = "Objective Difference",
                    breaks = c(1, 2, 3),
                    labels = c("< 2", "NA", "> 2")) +
  geom_point(aes(color = factor(Algorithm, levels = Names),
                 shape = factor(Convergence, levels = c(0, 1)),
                 size = factor(Obj_diff, levels = c(0, 1, 3)),
                 stroke = factor(Obj_diff, levels = c(0, 1, 3)),
                 fill = factor(Obj_diff, levels = c(0, 1, 3)))) +
  coord_cartesian(ylim = c(mife, mafe)) +
  labs(title = "Function evaluations for data fitting",
       subtitle = paste(
         "Dashed lines represent mean benchmark times",
         "Visible points have an objective difference greater than 1e-4",
         sep = ";\n"),
       x = "Method", y = "Number of function evaluations") +
  guides(shape = guide_legend(order = 1,
                              override.aes = list(size = Sizes[c(2, 3)])),
         fill = guide_legend(order = 2,
                             override.aes = list(size = Sizes[c(2, 3)],
                                                 shape = c(21, 21),
                                                 fill = Fills[c(2, 3)]))) +
  theme_bw() +
  theme(panel.border = element_blank(),
        plot.title = element_text(size = 23),
        plot.subtitle = element_text(size = 16,
                                     margin = margin(5, 5, 10, 5, "pt")),
        axis.text.x = element_text(size = 16, angle = 90,
                                   vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20,
                                    margin = margin(5, 10, 5, 5, "pt")),
        legend.position = "right",
        legend.box = "vertical",
        legend.direction = "vertical",
        legend.background = element_rect(fill = "transparent"),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 13))
ggsave("paper_analysis/image_files/bm_fit_fe.png", width = 16, height = 9)
