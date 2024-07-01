library("fddm")
library("microbenchmark")

ll_ft_SWSE_17 <- function(pars, rt, resp, truth, err_tol) {
  v <- numeric(length(rt))
  v[truth == "upper"] <- pars[[1]]
  v[truth == "lower"] <- pars[[2]]
  dens <- dfddm(rt = rt, response = resp, a = pars[[3]], v = v,
                t0 = pars[[4]], w = pars[[5]], sv = pars[[6]],
                err_tol = err_tol,
                log = TRUE, switch_mech = "eff_rt", switch_thresh = 0.8,
                n_terms_small = "SWSE", summation_small = "2017")
  return( ifelse(any(!is.finite(dens)), 1e6, -sum(dens)) )
}





rt_fit <- function(data, id_idx = NULL, rt_idx = NULL, response_idx = NULL,
                   truth_idx = NULL, response_upper = NULL, err_tol = 1e-6,
                   ctrl_list = list(eval.max = 1000, iter.max = 1000),
                   times = 100, unit = "ns") {

  # Format data for fitting
  if (all(is.null(id_idx), is.null(rt_idx), is.null(response_idx),
      is.null(truth_idx), is.null(response_upper))) {
    df <- data # assume input data is already formatted
  } else {
    if(any(data[,rt_idx] < 0)) {
      stop("Input data contains negative response times; fit will not be run.")
    }
    if(any(is.na(data[,response_idx]))) {
      stop("Input data contains invalid responses (NA); fit will not be run.")
    }

    nr <- nrow(data)
    df <- data.frame(id = character(nr),
                     rt = double(nr),
                     response = character(nr),
                     truth = character(nr),
                     stringsAsFactors = FALSE)

    if (!is.null(id_idx)) { # relabel identification tags
      for (i in seq_along(id_idx)) {
        idi <- unique(data[,id_idx[i]])
        for (j in seq_along(idi)) {
          df[["id"]][data[,id_idx[i]] == idi[j]] <- paste(
            df[["id"]][data[,id_idx[i]] == idi[j]], idi[j], sep = " ")
        }
      }
      df[["id"]] <- trimws(df[["id"]], which = "left")
    }

    df[["rt"]] <- as.double(data[,rt_idx])

    df[["response"]] <- "lower"
    df[["response"]][data[,response_idx] == response_upper] <- "upper"

    df[["truth"]] <- "lower"
    df[["truth"]][data[,truth_idx] == response_upper] <- "upper"
  }

  # Preliminaries
  ids <- unique(df[["id"]])
  nids <- max(length(ids), 1) # if inds is null, there is only one individual

  init_vals <- data.frame(v1 = c( 0,  10, -.5,  0,  0,  0,  0,  0,  0,   0,  0),
                          v0 = c( 0, -10,  .5,  0,  0,  0,  0,  0,  0,   0,  0),
                          a  = c( 1,   1,   1, .5,  5,  1,  1,  1,  1,   1,  1),
                          t0 = c( 0,   0,   0,  0,  0,  0,  0,  0,  0,   0,  0),
                          w  = c(.5,  .5,  .5, .5, .5, .5, .5, .2, .8,  .5, .5),
                          sv = c( 1,   1,   1,  1,  1,  1,  1,  1,  1, .05,  5))
  ninit_vals <- nrow(init_vals)
  ddm_inits <- data.frame(vint = c( 0, 10, -.5,  0,  0,  0,  0,  0,  0,   0,  0),
                          vdif = c( 0,  0,   0,  0,  0,  0,  0,  0,  0,   0,  0),
                          a    = c( 1,  1,   1, .5,  5,  1,  1,  1,  1,   1,  1),
                          t0   = c( 0,  0,   0,  0,  0,  0,  0,  0,  0,   0,  0),
                          w    = c(.5, .5,  .5, .5, .5, .5, .5, .2, .8,  .5, .5),
                          sv   = c( 1,  1,   1,  1,  1,  1,  1,  1,  1, .05,  5))

  algo_names <- c("manual_ll", "ddm")
  nalgos <- length(algo_names)
  ni <- nalgos*ninit_vals

  # Initilize the result dataframe
  cnames <- c("ID", "Algorithm", "Convergence", "Objective", "Iterations",
              "FuncEvals", "GradEvals", "BmTime")
  res <- data.frame(matrix(ncol = length(cnames),
                           nrow = nids * ninit_vals * nalgos))
  colnames(res) <- cnames

  # label the result dataframe
  res[["ID"]] <- rep(ids, each = ni) # label individuals
  res[["Algorithm"]] <- rep(algo_names, each = ninit_vals) # label algorithms

  # Loop through each individual
  for (i in 1:nids) {
    # extract data for id i
    dfi <- df[df[["id"]] == ids[i], ]
    rti <- dfi[["rt"]]
    respi <- dfi[["response"]]
    truthi <- dfi[["truth"]]

    # starting value for t0 must be smaller than the smallest rt
    min_rti <- min(rti)
    t0_lo <- 0.01*min_rti
    t0_me <- 0.50*min_rti
    t0_hi <- 0.99*min_rti
    init_vals[["t0"]] <- c(rep(t0_me, 5), t0_lo, t0_hi, rep(t0_me, 4))
    ddm_inits[["t0"]] <- c(rep(t0_me, 5), t0_lo, t0_hi, rep(t0_me, 4))

    # loop through all of the starting values
    for (j in 1:ninit_vals) {
      # get number of evaluations
      temp <- nlminb(init_vals[j, ], ll_ft_SWSE_17,
                     rt = rti, resp = respi, truth = truthi, err_tol = err_tol,
                     # limits:   vu,   vl,   a,      t0, w,  sv
                     lower = c(-Inf, -Inf, .01,       0, 0,   0),
                     upper = c( Inf,  Inf, Inf, min_rti, 1, Inf),
                     control = ctrl_list)
      res[["Convergence"]][(i-1)*ni+0*ninit_vals+j] <- temp[["convergence"]]
      res[["Objective"]][(i-1)*ni+0*ninit_vals+j] <- temp[["objective"]]
      res[["Iterations"]][(i-1)*ni+0*ninit_vals+j] <- temp[["iterations"]]
      res[["FuncEvals"]][(i-1)*ni+0*ninit_vals+j] <- temp[["evaluations"]][[1]]

      temp <- ddm(drift = rt + response ~ truth,
                  boundary = ~ 1, ndt = ~ 1, bias = ~ 1, sv = ~ 1, data = dfi,
                  args_optim = list(init = ddm_inits[j, ],
                                # limits:  vint, vdif,   a,      t0, w,  sv
                                lo_bds = c(-Inf, -Inf, .01,       0, 0,   0),
                                up_bds = c( Inf,  Inf, Inf, min_rti, 1, Inf),
                                control = ctrl_list))[["optim_info"]][["value"]]
      res[["Convergence"]][(i-1)*ni+1*ninit_vals+j] <- temp[["convergence"]]
      res[["Objective"]][(i-1)*ni+1*ninit_vals+j] <- temp[["objective"]]
      res[["Iterations"]][(i-1)*ni+1*ninit_vals+j] <- temp[["iterations"]]
      res[["FuncEvals"]][(i-1)*ni+1*ninit_vals+j] <- temp[["evaluations"]][[1]]
      res[["GradEvals"]][(i-1)*ni+1*ninit_vals+j] <- temp[["evaluations"]][[2]]

      # microbenchmark
      mbm <- microbenchmark(
        manual_ll = nlminb(init_vals[j, ], ll_ft_SWSE_17,
                           rt = rti, resp = respi, truth = truthi,
                           err_tol = err_tol,
                           # limits:   vu,   vl,   a,      t0, w,  sv
                           lower = c(-Inf, -Inf, .01,       0, 0,   0),
                           upper = c( Inf,  Inf, Inf, min_rti, 1, Inf),
                           control = ctrl_list),

        ddm = ddm(drift = rt + response ~ truth,
                  boundary = ~ 1, ndt = ~ 1, bias = ~ 1, sv = ~ 1, data = dfi,
                  args_optim = list(init = ddm_inits[j, ],
                                  # limits:  vint, vdif,   a,      t0, w,  sv
                                  lo_bds = c(-Inf, -Inf, .01,       0, 0,   0),
                                  up_bds = c( Inf,  Inf, Inf, min_rti, 1, Inf),
                                  control = ctrl_list)),
        times = times, unit = unit
      )
      for (k in 1:nalgos) {
        res[["BmTime"]][(i-1)*ni+(k-1)*ninit_vals+j] <- median(
          mbm[mbm[["expr"]] == algo_names[k], 2])
      }
    }
  }
  return(res)
}




###############################################################
data(med_dec, package = "fddm")
med_dec <- med_dec[which(med_dec[["rt"]] >= 0), ]
fit <- rt_fit(med_dec, id_idx = c(2,1), rt_idx = 8, response_idx = 7,
              truth_idx = 5, response_upper = "blast", err_tol = 1e-6,
              times = 5, unit = "ns")

save(fit, compress = "xz", compression_level = 9,
     file = "unused_stuff/bm_fit_ddm_vs_ll.Rds")


fit[fit[["Convergence"]] != 0,]



library("reshape2")
library("ggplot2")

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

obj_diff_label <- function(y, df, col_name, mult = 1.15, upper_limit = NULL) {
  if (is.null(upper_limit)) {
    upper_limit <- max(df[[as.character(col_name)]])
  }
  return(
    data.frame(
      y = mult * upper_limit,
      label = paste(sum(y > 0, na.rm = TRUE))
    )
  )
}

Names <- c("manual_ll", "ddm")
Color <- c("#92c639", "#ac8053")
Outline <- c("#92c639", "#ac8053")
Shape <- c(21, 25)
Sizes <- c(0, 3, 3)
Stroke <- c(0, 1, 1)
Fills <- c("#ffffff00", "#ffffff00", "#80808099")

# load data, will be in the variable 'fit'
load("/home/kendal/Documents/phd/henrik/fddm/unused_stuff/bm_fit_ddm_vs_ll.Rds")
fit <- fit_prep(fit)



fit_mbm <- melt(fit, id.vars = c("Algorithm", "Convergence", "Obj_diff"),
                measure.vars = "BmTime", value.name = "BmTime")[,-4]

ggplot(fit_mbm, aes(x = factor(Algorithm, levels = Names),
                    y = BmTime)) +
  geom_violin(trim = TRUE, alpha = 0.5,
              aes(color = factor(Algorithm, levels = Names),
                  fill = factor(Algorithm, levels = Names))) +
  geom_boxplot(width = 0.2, outlier.shape = NA,
               fill = "white", alpha = 0.4,
               aes(color = factor(Algorithm, levels = Names))) +
  stat_summary(fun = mean, geom = "errorbar",
               aes(ymax = after_stat(y), ymin = after_stat(y)),
               width = .5, linetype = "dashed",
               color = Color) +
  stat_summary(aes(y = Obj_diff, color = factor(Algorithm, levels = Names)),
               fun.data = obj_diff_label,
               fun.args = list(fit, "BmTime", 1.075, max(fit[["BmTime"]])),
               geom = "label",
               hjust = 0.5,
               vjust = 0.9) +
  scale_x_discrete(labels = c("manual_ll", "ddm")) +
  scale_color_manual(values = Outline, guide = "none") +
  scale_fill_manual(values = Color, guide = "none") +
  scale_shape_manual(values = Shape,
                     name = "Convergence Code",
                     breaks = c(0, 1),
                     labels = c("Success", "Failure")) +
  scale_size_manual(values = Sizes, guide = "none") +
  scale_discrete_manual(aesthetics = "stroke", values = Stroke, guide = "none") +
  ggnewscale::new_scale_fill() +
  scale_fill_manual(values = Fills,
                    name = paste("Difference in", "Log-likelihood", "from MLE",
                                 sep = "\n"),
                    breaks = c(1, 2, 3),
                    labels = c("< 2", "NA", "> 2")) +
  geom_point(aes(color = factor(Algorithm, levels = Names),
                 shape = factor(Convergence, levels = c(0, 1)),
                 size = factor(Obj_diff, levels = c(0, 1, 3)),
                 stroke = factor(Obj_diff, levels = c(0, 1, 3)),
                 fill = factor(Obj_diff, levels = c(0, 1, 3)))) +
  labs(x = "Implementation", y = "Time (milliseconds)") +
  guides(shape = guide_legend(order = 1,
                              override.aes = list(size = Sizes[c(2, 3)])),
         fill = guide_legend(order = 2,
                             override.aes = list(size = Sizes[c(2, 3)],
                                                 shape = c(21, 21),
                                                 fill = Fills[c(2, 3)]))) +
  theme_bw() +
  theme(panel.border = element_blank(),
        axis.text.x = element_text(size = 16, angle = 90,
                                   vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 20,
                                    margin = margin(15, 0, 0, 0)),
        axis.title.y = element_text(size = 20,
                                    margin = margin(0, 10, 0, 0)),
        legend.position = "right",
        legend.box = "vertical",
        legend.direction = "vertical",
        legend.background = element_rect(fill = "transparent"),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 13))






fit_fev <- melt(fit, id.vars = c("Algorithm", "Convergence", "Obj_diff"),
                measure.vars = "FuncEvals", value.name = "FuncEvals")[,-4]

ggplot(fit_fev, aes(x = factor(Algorithm, levels = Names),
                              y = FuncEvals)) +
  geom_violin(trim = TRUE, alpha = 0.5,
              aes(color = factor(Algorithm, levels = Names),
                  fill = factor(Algorithm, levels = Names))) +
  geom_boxplot(width = 0.2, outlier.shape = NA,
               fill = "white", alpha = 0.4,
               aes(color = factor(Algorithm, levels = Names))) +
  stat_summary(fun = mean, geom = "errorbar",
               aes(ymax = after_stat(y), ymin = after_stat(y)),
               width = .5, linetype = "dashed",
               color = Color) +
  stat_summary(aes(y = Obj_diff, color = factor(Algorithm, levels = Names)),
               fun.data = obj_diff_label,
               fun.args = list(fit, "FuncEvals", 1.075),
               geom = "label",
               hjust = 0.5,
               vjust = 0.9) +
  scale_x_discrete(labels = c("manual_ll", "ddm")) +
  scale_color_manual(values = Outline, guide = "none") +
  scale_fill_manual(values = Color, guide = "none") +
  scale_shape_manual(values = Shape,
                     name = "Convergence Code",
                     breaks = c(0, 1),
                     labels = c("Success", "Failure")) +
  scale_size_manual(values = Sizes, guide = "none") +
  scale_discrete_manual(aesthetics = "stroke", values = Stroke, guide = "none") +
  ggnewscale::new_scale_fill() +
  scale_fill_manual(values = Fills,
                    name = paste("Difference in", "Log-likelihood", "from MLE",
                                 sep = "\n"),
                    breaks = c(1, 2, 3),
                    labels = c("< 2", "NA", "> 2")) +
  geom_point(aes(color = factor(Algorithm, levels = Names),
                 shape = factor(Convergence, levels = c(0, 1)),
                 size = factor(Obj_diff, levels = c(0, 1, 3)),
                 stroke = factor(Obj_diff, levels = c(0, 1, 3)),
                 fill = factor(Obj_diff, levels = c(0, 1, 3)))) +
  labs(x = "Implementation", y = "Number of function evaluations") +
  guides(shape = guide_legend(order = 1,
                              override.aes = list(size = Sizes[c(2, 3)])),
         fill = guide_legend(order = 2,
                             override.aes = list(size = Sizes[c(2, 3)],
                                                 shape = c(21, 21),
                                                 fill = Fills[c(2, 3)]))) +
  theme_bw() +
  theme(panel.border = element_blank(),
        axis.text.x = element_text(size = 16, angle = 90,
                                   vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 20,
                                    margin = margin(15, 0, 0, 0)),
        axis.title.y = element_text(size = 20,
                                    margin = margin(0, 10, 0, 0)),
        legend.position = "right",
        legend.box = "vertical",
        legend.direction = "vertical",
        legend.background = element_rect(fill = "transparent"),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 13))

















############# gradient... doesn't work yet

onep <- med_dec[med_dec[["id"]] == "2" & med_dec[["group"]] == "experienced", ]
min_rti <- min(onep[["rt"]])
ctrl_list <- list(eval.max = 1000, iter.max = 1000)

onep[["resp"]] <- ifelse(onep[["response"]] == "blast", "upper", "lower")
onep[["truth"]] <- ifelse(onep[["classification"]] == "blast", "upper", "lower")




grad_ft_SWSE_17 <- function(pars, rt, resp, truth, err_tol) {
  v <- numeric(length(rt))
  v[truth == "upper"] <- pars[[1]]
  v[truth == "lower"] <- pars[[2]]
  dens <- dfddm(rt = rt, response = resp, a = pars[[3]], v = v,
                t0 = pars[[4]], w = pars[[5]], sv = pars[[6]],
                err_tol = err_tol,
                log = FALSE, switch_mech = "eff_rt", switch_thresh = 0.8,
                n_terms_small = "SWSE", summation_small = "2017")
  grad <- c(
    -1 * sum(dv_dfddm(rt = rt, response = resp, a = pars[[3]], v = v,
                      t0 = pars[[4]], w = pars[[5]], sv = pars[[6]],
                      err_tol = err_tol) / dens),
    sum(dv_dfddm(rt = rt, response = resp, a = pars[[3]], v = v, t0 = pars[[4]],
                 w = pars[[5]], sv = pars[[6]], err_tol = err_tol) / dens),
    -1 * sum(da_dfddm(rt = rt, response = resp, a = pars[[3]], v = v,
                      t0 = pars[[4]], w = pars[[5]], sv = pars[[6]],
                      err_tol = err_tol) / dens),
    -1 * sum(dt0_dfddm(rt = rt, response = resp, a = pars[[3]], v = v,
                       t0 = pars[[4]], w = pars[[5]], sv = pars[[6]],
                       err_tol = err_tol) / dens),
    -1 * sum(dw_dfddm(rt = rt, response = resp, a = pars[[3]], v = v,
                      t0 = pars[[4]], w = pars[[5]], sv = pars[[6]],
                      err_tol = err_tol) / dens),
    -1 * sum(dsv_dfddm(rt = rt, response = resp, a = pars[[3]], v = v,
                       t0 = pars[[4]], w = pars[[5]], sv = pars[[6]],
                       err_tol = err_tol) / dens)
  )
  return(grad)
}

grad_ft_SWSE_17(c(0, 0, 1, 0, 0.5, 0), onep[["rt"]], onep[["resp"]], onep[["truth"]], 1e-6)

nlminb(c(0, 0, 1, 0, 0.5, 0), ll_ft_SWSE_17,
       rt = onep[["rt"]], resp = onep[["resp"]], truth = onep[["truth"]],
       err_tol = 1e-6,
       # limits:   vu,   vl,   a,      t0, w,  sv
       lower = c(-Inf, -Inf, .01,       0, 0,   0),
       upper = c( Inf,  Inf, Inf, min_rti, 1, Inf),
       control = ctrl_list)
