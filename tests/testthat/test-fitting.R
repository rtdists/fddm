# context("Accuracy of fitting results across methods and established packages")
# library("rtdists")
# library("RWiener")
# source(system.file("extdata", "Kesselmeier_density.R", package = "fddm", mustWork = TRUE))
# 
# 
# ### See Known Errors (KE) at bottom for issues with Navarro's implementation
# ### Can read result dataframe in following line or rerun in lines 15-404
# fit <- readRDS(system.file("extdata", "test-fitting-fit.Rds", package = "fddm", mustWork = TRUE))
# 
# 
# ######################## Log-Likelihood functions ##############################
# ################################################################################
# 
# # ll_fs_Fos_17 <- function(pars, rt, resp, truth, t0, err_tol) {
# #   rt1 <- rt[truth == 1]
# #   rt0 <- rt[truth == 0]
# #   resp1 <- resp[truth == 1]
# #   resp0 <- resp[truth == 0]
# #   
# #   # the truth is "upper" so use v1
# #   dens1 <- dfddm(rt = rt1, response = resp1, a = pars[[3]],
# #                  v = pars[[1]], t0 = t0, w = pars[[4]],
# #                  log = FALSE, n_terms_small = "Foster",
# #                  summation_small = "2017", scale = "small",
# #                  err_tol = err_tol)
# #   # the truth is "lower" so use v0
# #   dens0 <- dfddm(rt = rt0, response = resp0, a = pars[[3]],
# #                  v = pars[[2]], t0 = t0, w = pars[[4]],
# #                  log = FALSE, n_terms_small = "Foster",
# #                  summation_small = "2017", scale = "small",
# #                  err_tol = err_tol)
# #   
# #   densities <- c(dens1, dens0)
# #   if (any(densities == 0)) return(1e6)
# #   return(-sum(log(densities)))
# # }
# # 
# # ll_fs_Kes_17 <- function(pars, rt, resp, truth, t0, err_tol) {
# #   rt1 <- rt[truth == 1]
# #   rt0 <- rt[truth == 0]
# #   resp1 <- resp[truth == 1]
# #   resp0 <- resp[truth == 0]
# #   
# #   # the truth is "upper" so use v1
# #   dens1 <- dfddm(rt = rt1, response = resp1, a = pars[[3]],
# #                  v = pars[[1]], t0 = t0, w = pars[[4]],
# #                  log = FALSE, n_terms_small = "Kesselmeier",
# #                  summation_small = "2017", scale = "small",
# #                  err_tol = err_tol)
# #   # the truth is "lower" so use v0
# #   dens0 <- dfddm(rt = rt0, response = resp0, a = pars[[3]],
# #                  v = pars[[2]], t0 = t0, w = pars[[4]],
# #                  log = FALSE, n_terms_small = "Kesselmeier",
# #                  summation_small = "2017", scale = "small",
# #                  err_tol = err_tol)
# #   
# #   densities <- c(dens1, dens0)
# #   if (any(densities == 0)) return(1e6)
# #   return(-sum(log(densities)))
# # }
# # 
# # ll_fs_Nav_17 <- function(pars, rt, resp, truth, t0, err_tol) {
# #   rt1 <- rt[truth == 1]
# #   rt0 <- rt[truth == 0]
# #   resp1 <- resp[truth == 1]
# #   resp0 <- resp[truth == 0]
# #   
# #   # the truth is "upper" so use v1
# #   dens1 <- dfddm(rt = rt1, response = resp1, a = pars[[3]],
# #                  v = pars[[1]], t0 = t0, w = pars[[4]],
# #                  log = FALSE, n_terms_small = "Navarro",
# #                  summation_small = "2017", scale = "small",
# #                  err_tol = err_tol)
# #   # the truth is "lower" so use v0
# #   dens0 <- dfddm(rt = rt0, response = resp0, a = pars[[3]],
# #                  v = pars[[2]], t0 = t0, w = pars[[4]],
# #                  log = FALSE, n_terms_small = "Navarro",
# #                  summation_small = "2017", scale = "small",
# #                  err_tol = err_tol)
# #   
# #   densities <- c(dens1, dens0)
# #   if (any(densities <= 0)) return(1e6) # Nav can give (small) negative densities
# #   return(-sum(log(densities)))
# # }
# # 
# # ll_fb_Kes_17 <- function(pars, rt, resp, truth, t0, err_tol) {
# #   rt1 <- rt[truth == 1]
# #   rt0 <- rt[truth == 0]
# #   resp1 <- resp[truth == 1]
# #   resp0 <- resp[truth == 0]
# #   
# #   # the truth is "upper" so use v1
# #   dens1 <- dfddm(rt = rt1, response = resp1, a = pars[[3]],
# #                  v = pars[[1]], t0 = t0, w = pars[[4]],
# #                  log = FALSE, n_terms_small = "Kesselmeier",
# #                  summation_small = "2017", scale = "both",
# #                  err_tol = err_tol)
# #   # the truth is "lower" so use v0
# #   dens0 <- dfddm(rt = rt0, response = resp0, a = pars[[3]],
# #                  v = pars[[2]], t0 = t0, w = pars[[4]],
# #                  log = FALSE, n_terms_small = "Kesselmeier",
# #                  summation_small = "2017", scale = "both",
# #                  err_tol = err_tol)
# #   
# #   densities <- c(dens1, dens0)
# #   if (any(densities == 0)) return(1e6)
# #   return(-sum(log(densities)))
# # }
# # 
# # ll_fb_Nav_17 <- function(pars, rt, resp, truth, t0, err_tol) {
# #   rt1 <- rt[truth == 1]
# #   rt0 <- rt[truth == 0]
# #   resp1 <- resp[truth == 1]
# #   resp0 <- resp[truth == 0]
# #   
# #   # the truth is "upper" so use v1
# #   dens1 <- dfddm(rt = rt1, response = resp1, a = pars[[3]],
# #                  v = pars[[1]], t0 = t0, w = pars[[4]],
# #                  log = FALSE, n_terms_small = "Navarro",
# #                  summation_small = "2017", scale = "both",
# #                  err_tol = err_tol)
# #   # the truth is "lower" so use v0
# #   dens0 <- dfddm(rt = rt0, response = resp0, a = pars[[3]],
# #                  v = pars[[2]], t0 = t0, w = pars[[4]],
# #                  log = FALSE, n_terms_small = "Navarro",
# #                  summation_small = "2017", scale = "both",
# #                  err_tol = err_tol)
# #   
# #   densities <- c(dens1, dens0)
# #   if (any(densities == 0)) return(1e6)
# #   return(-sum(log(densities)))
# # }
# # 
# # ll_RWiener <- function(pars, rt, resp, truth, t0) {
# #   rt1 <- rt[truth == 1]
# #   rt0 <- rt[truth == 0]
# #   resp1 <- resp[truth == 1]
# #   resp0 <- resp[truth == 0]
# #   
# #   # the truth is "upper" so use v1
# #   dens1 <- dwiener(rt1, resp = resp1, alpha = pars[[3]],
# #                    delta = pars[[1]], beta = pars[[4]], tau = t0)
# #   # the truth is "lower" so use v0
# #   dens0 <- dwiener(rt0, resp = resp0, alpha = pars[[3]],
# #                    delta = pars[[2]], beta = pars[[4]], tau = t0)
# #   
# #   densities <- c(dens1, dens0)
# #   if (any(densities == 0)) return(1e6)
# #   return(-sum(log(densities)))
# # }
# # 
# # ll_Kesselmeier <- function(pars, rt, resp, truth, t0, err_tol) {
# #   rt1 <- rt[truth == 1]
# #   rt0 <- rt[truth == 0]
# #   resp1 <- resp[truth == 1]
# #   resp0 <- resp[truth == 0]
# #   
# #   rt11 <- rt1[resp1 == 1]
# #   rt10 <- rt1[resp1 == 0]
# #   rt01 <- rt0[resp0 == 1]
# #   rt00 <- rt0[resp0 == 0]
# #   
# #   # truth is "upper" so use v1; resp is "upper" so use "upper" parameters
# #   dens11 <- fs14_R(t = rt11-t0, a = pars[[3]], v = -pars[[1]],
# #                    w = 1-pars[[4]], eps = err_tol)
# #   # truth is "upper" so use v1; resp is "lower" so use "lower" parameters
# #   dens10 <- fs14_R(t = rt10-t0, a = pars[[3]], v = pars[[1]],
# #                    w = pars[[4]], eps = err_tol)
# #   # truth is "lower" so use v0; resp is "upper" so use "upper" parameters
# #   dens01 <- fs14_R(t = rt01-t0, a = pars[[3]], v = -pars[[2]],
# #                    w = 1-pars[[4]], eps = err_tol)
# #   # truth is "lower" so use v0; resp is "lower" so use "lower" parameters
# #   dens00 <- fs14_R(t = rt00-t0, a = pars[[3]], v = pars[[2]],
# #                    w = pars[[4]], eps = err_tol)
# #   
# #   densities <- c(dens11, dens10, dens01, dens00)
# #   if (any(densities == 0)) return(1e6)
# #   return(-sum(log(densities)))
# # }
# # 
# # ll_RTDists <- function(pars, rt, resp, truth, t0) {
# #   rt1 <- rt[truth == 1]
# #   rt0 <- rt[truth == 0]
# #   resp1 <- resp[truth == 1]
# #   resp0 <- resp[truth == 0]
# #   
# #   # the truth is "upper" so use v1
# #   dens1 <- ddiffusion(rt1, resp1, a = pars[[3]],
# #                       v = pars[[1]], z = pars[[4]]*pars[[3]], t0 = t0)
# #   # the truth is "lower" so use v0
# #   dens0 <- ddiffusion(rt0, resp0, a = pars[[3]],
# #                       v = pars[[2]], z = pars[[4]]*pars[[3]], t0 = t0)
# #   
# #   densities <- c(dens1, dens0)
# #   if (any(densities == 0)) return(1e6)
# #   return(-sum(log(densities)))
# # }
# # 
# # 
# # 
# # 
# # ################################################################################
# # ############################ Fitting function ##################################
# # ################################################################################
# # 
# # rt_fit <- function(data, ind_idx = NULL, rt_idx = NULL, response_idx = NULL,
# #                    truth_idx = NULL, response_upper = NULL,
# #                    stvals = NULL, t0 = 1e-4, err_tol = 1e-6,
# #                    times = 1000, unit = "us")
# # {
# #   # Format data for fitting
# #   if (all(is.null(ind_idx), is.null(rt_idx), is.null(response_idx),
# #           is.null(truth_idx), is.null(response_upper))) {
# #     df <- data # assume input data is already formatted
# #   } else {
# #     nr <- nrow(data)
# #     df <- data.frame(ind = integer(nr),
# #                      rt = double(nr),
# #                      response = integer(nr),
# #                      rresponse = character(nr),
# #                      truth = integer(nr))
# #     
# #     df$ind <- as.integer(data[,ind_idx])
# #     df$rt <- as.double(data[,rt_idx])
# #     
# #     df$response <- as.integer(0)
# #     df$response[data[,response_idx] == response_upper] <- as.integer(1)
# #     df$rresponse <- "lower"
# #     df$rresponse[data[,response_idx] == response_upper] <- "upper"
# #     
# #     df$truth <- as.integer(0)
# #     df$truth[data[,truth_idx] == response_upper] <- as.integer(1)
# #   }
# #   
# #   # Preliminaries
# #   df <- subset(df, rt > t0) # remove any respnose times below t0
# #   inds <- sort(unique(df$ind))
# #   ninds <- max(length(inds), 1) # if inds is null, there is only one individual
# #   
# #   if (is.null(stvals)) {
# #     stvals <- list(#      v1,   v0,   a,   w
# #       list(   0,    0,   1, 0.5), # defaults
# #       list( 0.4, -0.4,   1, 0.5), # vary v1 & v0
# #       list(-0.4,  0.4,   1, 0.5), # vary v1 & v0
# #       list(   0,    0, 0.9, 0.5), # vary a; <.6 breaks Kesselmeier, <.9 breaks fs_Nav_17
# #       list(   0,    0, 2.5, 0.5), # vary a
# #       list(   0,    0,   1, 0.2), # vary w
# #       list(   0,    0,   1, 0.8)  # vary w
# #     )
# #   }
# #   nstvals <- length(stvals)
# #   stvals_mat <- matrix(unlist(stvals), ncol = 4, nrow = nstvals, byrow = TRUE)
# #   
# #   algo_names <- c(rep("fs_Fos_17", nstvals),
# #                   rep("fs_Kes_17", nstvals), rep("fs_Nav_17", nstvals),
# #                   rep("fb_Kes_17", nstvals), rep("fb_Nav_17", nstvals),
# #                   rep("RWiener", nstvals), rep("Kesselmeier", nstvals),
# #                   rep("rtdists", nstvals))
# #   nalgos <- length(unique(algo_names))
# #   ni <- nalgos*nstvals
# #   
# #   # Initilize the result dataframe
# #   res <- data.frame(matrix(ncol = 12, nrow = ninds*nstvals*nalgos))
# #   colnames(res) <- c("ind", "Algorithm", "AlgoCalls", "BmTime",
# #                      "init_v1", "init_v0", "init_a", "init_w",
# #                      "fit_v1", "fit_v0", "fit_a", "fit_w")
# #   
# #   # Fill in known values
# #   res[,2] <- algo_names # label algorithms
# #   res[,5] <- stvals_mat[,1] # label initial v1
# #   res[,6] <- stvals_mat[,2] # label initial v0
# #   res[,7] <- stvals_mat[,3] # label initial a
# #   res[,8] <- stvals_mat[,4] # label initial w
# #   
# #   # Loop through each individual and starting values
# #   for (i in 1:ninds) {
# #     dfi <- subset(df, ind == inds[i])
# #     res[((i-1)*ni+1):(i*ni), 1] <- inds[i] # label individuals
# #     rti <- dfi$rt
# #     respi <- dfi$response
# #     rrespi <- dfi$rresponse
# #     truthi <- dfi$truth
# #     for (j in 1:nstvals) {
# #       # Perform optimization and record results
# #       temp <- nlminb(stvals[[j]], ll_fs_Fos_17,
# #                      rt = rti, resp = respi, truth = truthi,
# #                      t0 = t0, err_tol = err_tol,
# #                      lower = c(-Inf, -Inf,   0, 0),
# #                      upper = c( Inf,  Inf, Inf, 1))
# #       res[(i-1)*ni+0*nstvals+j, 3] <- temp$evaluations[[1]]
# #       res[(i-1)*ni+0*nstvals+j, 9:12] <- temp$par
# #       
# #       temp <- nlminb(stvals[[j]], ll_fs_Kes_17,
# #                      rt = rti, resp = respi, truth = truthi,
# #                      t0 = t0, err_tol = err_tol,
# #                      lower = c(-Inf, -Inf,   0, 0),
# #                      upper = c( Inf,  Inf, Inf, 1))
# #       res[(i-1)*ni+1*nstvals+j, 3] <- temp$evaluations[[1]]
# #       res[(i-1)*ni+1*nstvals+j, 9:12] <- temp$par
# #       
# #       temp <- nlminb(stvals[[j]], ll_fs_Nav_17,
# #                      rt = rti, resp = respi, truth = truthi,
# #                      t0 = t0, err_tol = err_tol,
# #                      lower = c(-Inf, -Inf,   0, 0),
# #                      upper = c( Inf,  Inf, Inf, 1))
# #       res[(i-1)*ni+2*nstvals+j, 3] <- temp$evaluations[[1]]
# #       res[(i-1)*ni+2*nstvals+j, 9:12] <- temp$par
# #       
# #       temp <- nlminb(stvals[[j]], ll_fb_Kes_17,
# #                      rt = rti, resp = respi, truth = truthi,
# #                      t0 = t0, err_tol = err_tol,
# #                      lower = c(-Inf, -Inf,   0, 0),
# #                      upper = c( Inf,  Inf, Inf, 1))
# #       res[(i-1)*ni+3*nstvals+j, 3] <- temp$evaluations[[1]]
# #       res[(i-1)*ni+3*nstvals+j, 9:12] <- temp$par
# #       
# #       temp <- nlminb(stvals[[j]], ll_fb_Nav_17,
# #                      rt = rti, resp = respi, truth = truthi,
# #                      t0 = t0, err_tol = err_tol,
# #                      lower = c(-Inf, -Inf,   0, 0),
# #                      upper = c( Inf,  Inf, Inf, 1))
# #       res[(i-1)*ni+4*nstvals+j, 3] <- temp$evaluations[[1]]
# #       res[(i-1)*ni+4*nstvals+j, 9:12] <- temp$par
# #       
# #       temp <- nlminb(stvals[[j]], ll_RWiener,
# #                      rt = rti, resp = rrespi, truth = truthi,
# #                      t0 = t0,
# #                      lower = c(-Inf, -Inf,   0, 0),
# #                      upper = c( Inf,  Inf, Inf, 1))
# #       res[(i-1)*ni+5*nstvals+j, 3] <- temp$evaluations[[1]]
# #       res[(i-1)*ni+5*nstvals+j, 9:12] <- temp$par
# #       
# #       temp <- nlminb(stvals[[j]], ll_Kesselmeier,
# #                      rt = rti, resp = respi, truth = truthi,
# #                      t0 = t0, err_tol = err_tol,
# #                      lower = c(-Inf, -Inf,   0.01, 0), # lower bound for a > 0
# #                      upper = c( Inf,  Inf, Inf, 1))
# #       res[(i-1)*ni+6*nstvals+j, 3] <- temp$evaluations[[1]]
# #       res[(i-1)*ni+6*nstvals+j, 9:12] <- temp$par
# #       
# #       temp <- nlminb(stvals[[j]], ll_RTDists,
# #                      rt = rti, resp = rrespi, truth = truthi,
# #                      t0 = t0,
# #                      lower = c(-Inf, -Inf,   0, 0),
# #                      upper = c( Inf,  Inf, Inf, 1))
# #       res[(i-1)*ni+7*nstvals+j, 3] <- temp$evaluations[[1]]
# #       res[(i-1)*ni+7*nstvals+j, 9:12] <- temp$par
# #       
# #       # Benchmark optimizations
# #       mbm <- microbenchmark(
# #         fs_Fos_17 = nlminb(stvals[[j]], ll_fs_Fos_17,
# #                            rt = rti, resp = respi, truth = truthi,
# #                            t0 = t0, err_tol = err_tol,
# #                            lower = c(-Inf, -Inf,   0, 0),
# #                            upper = c( Inf,  Inf, Inf, 1)),
# #         fs_Kes_17 = nlminb(stvals[[j]], ll_fs_Kes_17,
# #                            rt = rti, resp = respi, truth = truthi,
# #                            t0 = t0, err_tol = err_tol,
# #                            lower = c(-Inf, -Inf,   0, 0),
# #                            upper = c( Inf,  Inf, Inf, 1)),
# #         fs_Nav_17 = nlminb(stvals[[j]], ll_fs_Nav_17,
# #                            rt = rti, resp = respi, truth = truthi,
# #                            t0 = t0, err_tol = err_tol,
# #                            lower = c(-Inf, -Inf,   0, 0),
# #                            upper = c( Inf,  Inf, Inf, 1)),
# #         fb_Kes_17 = nlminb(stvals[[j]], ll_fb_Kes_17,
# #                            rt = rti, resp = respi, truth = truthi,
# #                            t0 = t0, err_tol = err_tol,
# #                            lower = c(-Inf, -Inf,   0, 0),
# #                            upper = c( Inf,  Inf, Inf, 1)),
# #         fb_Nav_17 = nlminb(stvals[[j]], ll_fb_Nav_17,
# #                            rt = rti, resp = respi, truth = truthi,
# #                            t0 = t0, err_tol = err_tol,
# #                            lower = c(-Inf, -Inf,   0, 0),
# #                            upper = c( Inf,  Inf, Inf, 1)),
# #         RWiener = nlminb(stvals[[j]], ll_RWiener,
# #                          rt = rti, resp = rrespi, truth = truthi,
# #                          t0 = t0,
# #                          lower = c(-Inf, -Inf,   0, 0),
# #                          upper = c( Inf,  Inf, Inf, 1)),
# #         Kesselmeier = nlminb(stvals[[j]], ll_Kesselmeier,
# #                              rt = rti, resp = respi, truth = truthi,
# #                              t0 = t0, err_tol = err_tol,
# #                              lower = c(-Inf, -Inf,   0.01, 0), # lower bound for a > 0
# #                              upper = c( Inf,  Inf, Inf, 1)),
# #         rtdists = nlminb(stvals[[j]], ll_RTDists,
# #                          rt = rti, resp = rrespi, truth = truthi,
# #                          t0 = t0,
# #                          lower = c(-Inf, -Inf,   0, 0),
# #                          upper = c( Inf,  Inf, Inf, 1)),
# #         times = times, unit = unit
# #       )
# #       for (k in 1:nalgos) {
# #         res[(i-1)*ni+(k-1)*nstvals+j, 4] <- median(subset(mbm, expr == algo_names[k])$time)
# #       }
# #     }
# #   }
# #   return(res)
# # }
# # # Import raw data and rerun all fits
# # data(med_dec, package = "fddm")
# # fit <- rt_fit(med_dec, ind_idx = 1, rt_idx = 8, response_idx = 7, truth_idx = 5,
# #               response_upper = "blast", times = 10000, unit = "us")
# 
# 
# 
# 
# ### Check fits are all approximately equal
# # Define error tolerance
# eps = 1e-4
# 
# # Calculate differences and means
# fit$v1_dif <- rep(0, nrow(fit))
# fit$v0_dif <- rep(0, nrow(fit))
# fit$a_dif <- rep(0, nrow(fit))
# fit$w_dif <- rep(0, nrow(fit))
# 
# inds <- unique(fit$ind)
# ninds <- length(inds)
# algos <- unique(fit$Algorithm)
# nalgos <- length(algos)
# 
# avg_df <- data.frame(matrix(ncol = 6, nrow = ninds*nalgos))
# colnames(avg_df) <- c('ind', 'Algorithm', 'v1_avg', 'v0_avg', 'a_avg', 'w_avg')
# avg_df$Algorithm <- rep(algos, ninds)
# 
# start <- 1
# stop <- nrow(subset(fit, ind == inds[1] & Algorithm == algos[1]))
# for (i in 1:ninds) {
#   avg_df$ind[((i-1)*nalgos+1):(i*nalgos)] <- rep(inds[i], nalgos) # label ind
#   for (j in 1:nalgos) {
#     tmp <- subset(fit, ind == inds[i] & Algorithm == algos[j])[,8:11]
#     avg_df[(i-1)*nalgos+j,3:6] <- colMeans(tmp) # collect means
#     for (k in 1:(stop-1)) {
#       fit[start+k, 12:15] <- abs(fit[start, 8:11] - fit[(start+k), 8:11])
#     }
#     start = start + stop
#   }
# }
# 
# 
# # Subset
# Foster <- subset(fit, Algorithm == "fs_Fos_17")
# Kesselmeier_s <- subset(fit, Algorithm == "fs_Kes_17")
# Navarro_s <- subset(fit, Algorithm == "fs_Nav_17")
# Kesselmeier_b <- subset(fit, Algorithm == "fb_Kes_17")
# Navarro_b <- subset(fit, Algorithm == "fb_Nav_17")
# RWiener <- subset(fit, Algorithm == "RWiener")
# Kesselmeier_R <- subset(fit, Algorithm == "Kesselmeier")
# rtdists <- subset(fit, Algorithm == "rtdists")
# 
# Foster_avg <- subset(avg_df, Algorithm == "fs_Fos_17")
# Kesselmeier_s_avg <- subset(avg_df, Algorithm == "fs_Kes_17")
# Navarro_s_avg <- subset(avg_df, Algorithm == "fs_Nav_17")
# Kesselmeier_b_avg <- subset(avg_df, Algorithm == "fb_Kes_17")
# Navarro_b_avg <- subset(avg_df, Algorithm == "fb_Nav_17")
# RWiener_avg <- subset(avg_df, Algorithm == "RWiener")
# Kesselmeier_R_avg <- subset(avg_df, Algorithm == "Kesselmeier")
# rtdists_avg <- subset(avg_df, Algorithm == "rtdists")
# 
# 
# # Test accuracy for fitted parameters: v_upper, v_lower, a, w
# test_that("Consistency within each method, using different starting values", {
#   expect_true(all(Foster[,13:16] < eps))
#   expect_true(all(Kesselmeier_s[,13:16] < eps))
#   expect_true(all(Navarro_s[,13:16] < eps))
#   expect_true(all(Kesselmeier_b[,13:16] < eps))
#   expect_true(all(Navarro_b[,13:16] < eps))
#   expect_true(all(RWiener[,13:16] < eps))
#   expect_true(all(Kesselmeier_R[,13:16] < eps))
#   expect_true(all(rtdists[,13:16] < eps))
# })
# 
# test_that("Accuracy across different methods", {
#   expect_true(all(abs(Foster_avg[,-(1:2)] - Kesselmeier_s_avg[,-(1:2)]) < eps))
#   expect_true(all(abs(Foster_avg[,-(1:2)] - Navarro_s_avg[,-(1:2)]) < eps))
#   expect_true(all(abs(Foster_avg[,-(1:2)] - Kesselmeier_b_avg[,-(1:2)]) < eps))
#   expect_true(all(abs(Foster_avg[,-(1:2)] - Navarro_b_avg[,-(1:2)]) < eps))
#   expect_true(all(abs(Foster_avg[,-(1:2)] - RWiener_avg[,-(1:2)]) < eps))
#   expect_true(all(abs(Foster_avg[,-(1:2)] - Kesselmeier_R_avg[,-(1:2)]) < eps))
#   expect_true(all(abs(Foster_avg[,-(1:2)] - rtdists_avg[,-(1:2)]) < eps))
# })
