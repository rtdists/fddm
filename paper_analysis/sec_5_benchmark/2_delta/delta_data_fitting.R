# Figure 5

# This file produces the results (and plots) that pertain to data fitting using
# the implementations that combine the SWSE small-time method with the
# Navarro et al. large-time method, as shown in
# Section 5.2 of the fddm paper:
# "Determining the Default Behavior of our Heuristic Switching Mechanism, delta"

library("fddm")
library("microbenchmark")
library("reshape2")
library("ggplot2")
save_dir <- "paper_analysis/sec_5_benchmark/2_delta/results/"
img_dir <- "paper_analysis/images/sec_5_benchmark/"



##### Log-Likelihood Functions #################################################
ll_fb_SWSE_17_0 <- function(pars, rt, resp, truth, err_tol) {
  v <- numeric(length(rt))
  v[truth == "upper"] <- pars[[1]]
  v[truth == "lower"] <- pars[[2]]
  dens <- dfddm(rt = rt, response = resp, a = pars[[3]], v = v,
                t0 = pars[[4]], w = pars[[5]], sv = pars[[6]], log = TRUE,
                n_terms_small = "SWSE", summation_small = "2017",
                scale = "both", max_terms_large = 0, err_tol = err_tol)
  return( ifelse(any(!is.finite(dens)), 1e6, -sum(dens)) )
}

ll_fb_SWSE_14_0 <- function(pars, rt, resp, truth, err_tol) {
  v <- numeric(length(rt))
  v[truth == "upper"] <- pars[[1]]
  v[truth == "lower"] <- pars[[2]]
  dens <- dfddm(rt = rt, response = resp, a = pars[[3]], v = v,
                t0 = pars[[4]], w = pars[[5]], sv = pars[[6]], log = TRUE,
                n_terms_small = "SWSE", summation_small = "2014",
                scale = "both", max_terms_large = 0, err_tol = err_tol)
  return( ifelse(any(!is.finite(dens)), 1e6, -sum(dens)) )
}

ll_fb_SWSE_17_1 <- function(pars, rt, resp, truth, err_tol) {
  v <- numeric(length(rt))
  v[truth == "upper"] <- pars[[1]]
  v[truth == "lower"] <- pars[[2]]
  dens <- dfddm(rt = rt, response = resp, a = pars[[3]], v = v,
                t0 = pars[[4]], w = pars[[5]], sv = pars[[6]], log = TRUE,
                n_terms_small = "SWSE", summation_small = "2017",
                scale = "both", max_terms_large = 1, err_tol = err_tol)
  return( ifelse(any(!is.finite(dens)), 1e6, -sum(dens)) )
}

ll_fb_SWSE_14_1 <- function(pars, rt, resp, truth, err_tol) {
  v <- numeric(length(rt))
  v[truth == "upper"] <- pars[[1]]
  v[truth == "lower"] <- pars[[2]]
  dens <- dfddm(rt = rt, response = resp, a = pars[[3]], v = v,
                t0 = pars[[4]], w = pars[[5]], sv = pars[[6]], log = TRUE,
                n_terms_small = "SWSE", summation_small = "2014",
                scale = "both", max_terms_large = 1, err_tol = err_tol)
  return( ifelse(any(!is.finite(dens)), 1e6, -sum(dens)) )
}

ll_fb_SWSE_17_2 <- function(pars, rt, resp, truth, err_tol) {
  v <- numeric(length(rt))
  v[truth == "upper"] <- pars[[1]]
  v[truth == "lower"] <- pars[[2]]
  dens <- dfddm(rt = rt, response = resp, a = pars[[3]], v = v,
                t0 = pars[[4]], w = pars[[5]], sv = pars[[6]], log = TRUE,
                n_terms_small = "SWSE", summation_small = "2017",
                scale = "both", max_terms_large = 2, err_tol = err_tol)
  return( ifelse(any(!is.finite(dens)), 1e6, -sum(dens)) )
}

ll_fb_SWSE_14_2 <- function(pars, rt, resp, truth, err_tol) {
  v <- numeric(length(rt))
  v[truth == "upper"] <- pars[[1]]
  v[truth == "lower"] <- pars[[2]]
  dens <- dfddm(rt = rt, response = resp, a = pars[[3]], v = v,
                t0 = pars[[4]], w = pars[[5]], sv = pars[[6]], log = TRUE,
                n_terms_small = "SWSE", summation_small = "2014",
                scale = "both", max_terms_large = 2, err_tol = err_tol)
  return( ifelse(any(!is.finite(dens)), 1e6, -sum(dens)) )
}

ll_fb_SWSE_17_3 <- function(pars, rt, resp, truth, err_tol) {
  v <- numeric(length(rt))
  v[truth == "upper"] <- pars[[1]]
  v[truth == "lower"] <- pars[[2]]
  dens <- dfddm(rt = rt, response = resp, a = pars[[3]], v = v,
                t0 = pars[[4]], w = pars[[5]], sv = pars[[6]], log = TRUE,
                n_terms_small = "SWSE", summation_small = "2017",
                scale = "both", max_terms_large = 3, err_tol = err_tol)
  return( ifelse(any(!is.finite(dens)), 1e6, -sum(dens)) )
}

ll_fb_SWSE_14_3 <- function(pars, rt, resp, truth, err_tol) {
  v <- numeric(length(rt))
  v[truth == "upper"] <- pars[[1]]
  v[truth == "lower"] <- pars[[2]]
  dens <- dfddm(rt = rt, response = resp, a = pars[[3]], v = v,
                t0 = pars[[4]], w = pars[[5]], sv = pars[[6]], log = TRUE,
                n_terms_small = "SWSE", summation_small = "2014",
                scale = "both", max_terms_large = 3, err_tol = err_tol)
  return( ifelse(any(!is.finite(dens)), 1e6, -sum(dens)) )
}

ll_fb_SWSE_17_4 <- function(pars, rt, resp, truth, err_tol) {
  v <- numeric(length(rt))
  v[truth == "upper"] <- pars[[1]]
  v[truth == "lower"] <- pars[[2]]
  dens <- dfddm(rt = rt, response = resp, a = pars[[3]], v = v,
                t0 = pars[[4]], w = pars[[5]], sv = pars[[6]], log = TRUE,
                n_terms_small = "SWSE", summation_small = "2017",
                scale = "both", max_terms_large = 4, err_tol = err_tol)
  return( ifelse(any(!is.finite(dens)), 1e6, -sum(dens)) )
}

ll_fb_SWSE_14_4 <- function(pars, rt, resp, truth, err_tol) {
  v <- numeric(length(rt))
  v[truth == "upper"] <- pars[[1]]
  v[truth == "lower"] <- pars[[2]]
  dens <- dfddm(rt = rt, response = resp, a = pars[[3]], v = v,
                t0 = pars[[4]], w = pars[[5]], sv = pars[[6]], log = TRUE,
                n_terms_small = "SWSE", summation_small = "2014",
                scale = "both", max_terms_large = 4, err_tol = err_tol)
  return( ifelse(any(!is.finite(dens)), 1e6, -sum(dens)) )
}

ll_fb_SWSE_17_5 <- function(pars, rt, resp, truth, err_tol) {
  v <- numeric(length(rt))
  v[truth == "upper"] <- pars[[1]]
  v[truth == "lower"] <- pars[[2]]
  dens <- dfddm(rt = rt, response = resp, a = pars[[3]], v = v,
                t0 = pars[[4]], w = pars[[5]], sv = pars[[6]], log = TRUE,
                n_terms_small = "SWSE", summation_small = "2017",
                scale = "both", max_terms_large = 5, err_tol = err_tol)
  return( ifelse(any(!is.finite(dens)), 1e6, -sum(dens)) )
}

ll_fb_SWSE_14_5 <- function(pars, rt, resp, truth, err_tol) {
  v <- numeric(length(rt))
  v[truth == "upper"] <- pars[[1]]
  v[truth == "lower"] <- pars[[2]]
  dens <- dfddm(rt = rt, response = resp, a = pars[[3]], v = v,
                t0 = pars[[4]], w = pars[[5]], sv = pars[[6]], log = TRUE,
                n_terms_small = "SWSE", summation_small = "2014",
                scale = "both", max_terms_large = 5, err_tol = err_tol)
  return( ifelse(any(!is.finite(dens)), 1e6, -sum(dens)) )
}

ll_fb_SWSE_17_6 <- function(pars, rt, resp, truth, err_tol) {
  v <- numeric(length(rt))
  v[truth == "upper"] <- pars[[1]]
  v[truth == "lower"] <- pars[[2]]
  dens <- dfddm(rt = rt, response = resp, a = pars[[3]], v = v,
                t0 = pars[[4]], w = pars[[5]], sv = pars[[6]], log = TRUE,
                n_terms_small = "SWSE", summation_small = "2017",
                scale = "both", max_terms_large = 6, err_tol = err_tol)
  return( ifelse(any(!is.finite(dens)), 1e6, -sum(dens)) )
}

ll_fb_SWSE_14_6 <- function(pars, rt, resp, truth, err_tol) {
  v <- numeric(length(rt))
  v[truth == "upper"] <- pars[[1]]
  v[truth == "lower"] <- pars[[2]]
  dens <- dfddm(rt = rt, response = resp, a = pars[[3]], v = v,
                t0 = pars[[4]], w = pars[[5]], sv = pars[[6]], log = TRUE,
                n_terms_small = "SWSE", summation_small = "2014",
                scale = "both", max_terms_large = 6, err_tol = err_tol)
  return( ifelse(any(!is.finite(dens)), 1e6, -sum(dens)) )
}

ll_fb_SWSE_17_7 <- function(pars, rt, resp, truth, err_tol) {
  v <- numeric(length(rt))
  v[truth == "upper"] <- pars[[1]]
  v[truth == "lower"] <- pars[[2]]
  dens <- dfddm(rt = rt, response = resp, a = pars[[3]], v = v,
                t0 = pars[[4]], w = pars[[5]], sv = pars[[6]], log = TRUE,
                n_terms_small = "SWSE", summation_small = "2017",
                scale = "both", max_terms_large = 7, err_tol = err_tol)
  return( ifelse(any(!is.finite(dens)), 1e6, -sum(dens)) )
}

ll_fb_SWSE_14_7 <- function(pars, rt, resp, truth, err_tol) {
  v <- numeric(length(rt))
  v[truth == "upper"] <- pars[[1]]
  v[truth == "lower"] <- pars[[2]]
  dens <- dfddm(rt = rt, response = resp, a = pars[[3]], v = v,
                t0 = pars[[4]], w = pars[[5]], sv = pars[[6]], log = TRUE,
                n_terms_small = "SWSE", summation_small = "2014",
                scale = "both", max_terms_large = 7, err_tol = err_tol)
  return( ifelse(any(!is.finite(dens)), 1e6, -sum(dens)) )
}


##### Benchmark Functions ######################################################
rt_fit <- function(data, id_idx = NULL, rt_idx = NULL, response_idx = NULL,
                   truth_idx = NULL, response_upper = NULL, err_tol = 1e-6,
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
      for (i in 1:length(id_idx)) {
        idi <- unique(data[,id_idx[i]])
        for (j in 1:length(idi)) {
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

  algo_names <- c("fb_SWSE_17_0", "fb_SWSE_14_0", "fb_SWSE_17_1", "fb_SWSE_14_1",
                  "fb_SWSE_17_2", "fb_SWSE_14_2", "fb_SWSE_17_3", "fb_SWSE_14_3",
                  "fb_SWSE_17_4", "fb_SWSE_14_4", "fb_SWSE_17_5", "fb_SWSE_14_5",
                  "fb_SWSE_17_6", "fb_SWSE_14_6", "fb_SWSE_17_7", "fb_SWSE_14_7")
  nalgos <- length(algo_names)
  ni <- nalgos*ninit_vals

  # Initilize the result dataframe
  cnames <- c("ID", "Algorithm", "Convergence", "Objective", "Iterations",
              "FuncEvals", "BmTime")
  res <- data.frame(matrix(ncol = length(cnames), nrow = nids*ninit_vals*nalgos))
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

    # loop through all of the starting values
    for (j in 1:ninit_vals) {
      # get number of evaluations
      temp <- nlminb(init_vals[j, ], ll_fb_SWSE_17_0,
                     rt = rti, resp = respi, truth = truthi, err_tol = err_tol,
                     # limits:   vu,   vl,   a,      t0, w,  sv
                     lower = c(-Inf, -Inf, .01,       0, 0,   0),
                     upper = c( Inf,  Inf, Inf, min_rti, 1, Inf))
      res[["Convergence"]][(i-1)*ni+0*ninit_vals+j] <- temp[["convergence"]]
      res[["Objective"]][(i-1)*ni+0*ninit_vals+j] <- temp[["objective"]]
      res[["Iterations"]][(i-1)*ni+0*ninit_vals+j] <- temp[["iterations"]]
      res[["FuncEvals"]][(i-1)*ni+0*ninit_vals+j] <- temp[["evaluations"]][[1]]

      temp <- nlminb(init_vals[j, ], ll_fb_SWSE_14_0,
                     rt = rti, resp = respi, truth = truthi, err_tol = err_tol,
                     # limits:   vu,   vl,   a,      t0, w,  sv
                     lower = c(-Inf, -Inf, .01,       0, 0,   0),
                     upper = c( Inf,  Inf, Inf, min_rti, 1, Inf))
      res[["Convergence"]][(i-1)*ni+1*ninit_vals+j] <- temp[["convergence"]]
      res[["Objective"]][(i-1)*ni+1*ninit_vals+j] <- temp[["objective"]]
      res[["Iterations"]][(i-1)*ni+1*ninit_vals+j] <- temp[["iterations"]]
      res[["FuncEvals"]][(i-1)*ni+1*ninit_vals+j] <- temp[["evaluations"]][[1]]

      temp <- nlminb(init_vals[j, ], ll_fb_SWSE_17_1,
                     rt = rti, resp = respi, truth = truthi, err_tol = err_tol,
                     # limits:   vu,   vl,   a,      t0, w,  sv
                     lower = c(-Inf, -Inf, .01,       0, 0,   0),
                     upper = c( Inf,  Inf, Inf, min_rti, 1, Inf))
      res[["Convergence"]][(i-1)*ni+2*ninit_vals+j] <- temp[["convergence"]]
      res[["Objective"]][(i-1)*ni+2*ninit_vals+j] <- temp[["objective"]]
      res[["Iterations"]][(i-1)*ni+2*ninit_vals+j] <- temp[["iterations"]]
      res[["FuncEvals"]][(i-1)*ni+2*ninit_vals+j] <- temp[["evaluations"]][[1]]

      temp <- nlminb(init_vals[j, ], ll_fb_SWSE_14_1,
                     rt = rti, resp = respi, truth = truthi, err_tol = err_tol,
                     # limits:   vu,   vl,   a,      t0, w,  sv
                     lower = c(-Inf, -Inf, .01,       0, 0,   0),
                     upper = c( Inf,  Inf, Inf, min_rti, 1, Inf))
      res[["Convergence"]][(i-1)*ni+3*ninit_vals+j] <- temp[["convergence"]]
      res[["Objective"]][(i-1)*ni+3*ninit_vals+j] <- temp[["objective"]]
      res[["Iterations"]][(i-1)*ni+3*ninit_vals+j] <- temp[["iterations"]]
      res[["FuncEvals"]][(i-1)*ni+3*ninit_vals+j] <- temp[["evaluations"]][[1]]

      temp <- nlminb(init_vals[j, ], ll_fb_SWSE_17_2,
                     rt = rti, resp = respi, truth = truthi, err_tol = err_tol,
                     # limits:   vu,   vl,   a,      t0, w,  sv
                     lower = c(-Inf, -Inf, .01,       0, 0,   0),
                     upper = c( Inf,  Inf, Inf, min_rti, 1, Inf))
      res[["Convergence"]][(i-1)*ni+4*ninit_vals+j] <- temp[["convergence"]]
      res[["Objective"]][(i-1)*ni+4*ninit_vals+j] <- temp[["objective"]]
      res[["Iterations"]][(i-1)*ni+4*ninit_vals+j] <- temp[["iterations"]]
      res[["FuncEvals"]][(i-1)*ni+4*ninit_vals+j] <- temp[["evaluations"]][[1]]

      temp <- nlminb(init_vals[j, ], ll_fb_SWSE_14_2,
                     rt = rti, resp = respi, truth = truthi, err_tol = err_tol,
                     # limits:   vu,   vl,   a,      t0, w,  sv
                     lower = c(-Inf, -Inf, .01,       0, 0,   0),
                     upper = c( Inf,  Inf, Inf, min_rti, 1, Inf))
      res[["Convergence"]][(i-1)*ni+5*ninit_vals+j] <- temp[["convergence"]]
      res[["Objective"]][(i-1)*ni+5*ninit_vals+j] <- temp[["objective"]]
      res[["Iterations"]][(i-1)*ni+5*ninit_vals+j] <- temp[["iterations"]]
      res[["FuncEvals"]][(i-1)*ni+5*ninit_vals+j] <- temp[["evaluations"]][[1]]

      temp <- nlminb(init_vals[j, ], ll_fb_SWSE_17_3,
                     rt = rti, resp = respi, truth = truthi, err_tol = err_tol,
                     # limits:   vu,   vl,   a,      t0, w,  sv
                     lower = c(-Inf, -Inf, .01,       0, 0,   0),
                     upper = c( Inf,  Inf, Inf, min_rti, 1, Inf))
      res[["Convergence"]][(i-1)*ni+6*ninit_vals+j] <- temp[["convergence"]]
      res[["Objective"]][(i-1)*ni+6*ninit_vals+j] <- temp[["objective"]]
      res[["Iterations"]][(i-1)*ni+6*ninit_vals+j] <- temp[["iterations"]]
      res[["FuncEvals"]][(i-1)*ni+6*ninit_vals+j] <- temp[["evaluations"]][[1]]

      temp <- nlminb(init_vals[j, ], ll_fb_SWSE_14_3,
                     rt = rti, resp = respi, truth = truthi, err_tol = err_tol,
                     # limits:   vu,   vl,   a,      t0, w,  sv
                     lower = c(-Inf, -Inf, .01,       0, 0,   0),
                     upper = c( Inf,  Inf, Inf, min_rti, 1, Inf))
      res[["Convergence"]][(i-1)*ni+7*ninit_vals+j] <- temp[["convergence"]]
      res[["Objective"]][(i-1)*ni+7*ninit_vals+j] <- temp[["objective"]]
      res[["Iterations"]][(i-1)*ni+7*ninit_vals+j] <- temp[["iterations"]]
      res[["FuncEvals"]][(i-1)*ni+7*ninit_vals+j] <- temp[["evaluations"]][[1]]

      temp <- nlminb(init_vals[j, ], ll_fb_SWSE_17_4,
                     rt = rti, resp = respi, truth = truthi, err_tol = err_tol,
                     # limits:   vu,   vl,   a,      t0, w,  sv
                     lower = c(-Inf, -Inf, .01,       0, 0,   0),
                     upper = c( Inf,  Inf, Inf, min_rti, 1, Inf))
      res[["Convergence"]][(i-1)*ni+8*ninit_vals+j] <- temp[["convergence"]]
      res[["Objective"]][(i-1)*ni+8*ninit_vals+j] <- temp[["objective"]]
      res[["Iterations"]][(i-1)*ni+8*ninit_vals+j] <- temp[["iterations"]]
      res[["FuncEvals"]][(i-1)*ni+8*ninit_vals+j] <- temp[["evaluations"]][[1]]

      temp <- nlminb(init_vals[j, ], ll_fb_SWSE_14_4,
                     rt = rti, resp = respi, truth = truthi, err_tol = err_tol,
                     # limits:   vu,   vl,   a,      t0, w,  sv
                     lower = c(-Inf, -Inf, .01,       0, 0,   0),
                     upper = c( Inf,  Inf, Inf, min_rti, 1, Inf))
      res[["Convergence"]][(i-1)*ni+9*ninit_vals+j] <- temp[["convergence"]]
      res[["Objective"]][(i-1)*ni+9*ninit_vals+j] <- temp[["objective"]]
      res[["Iterations"]][(i-1)*ni+9*ninit_vals+j] <- temp[["iterations"]]
      res[["FuncEvals"]][(i-1)*ni+9*ninit_vals+j] <- temp[["evaluations"]][[1]]

      temp <- nlminb(init_vals[j, ], ll_fb_SWSE_17_5,
                     rt = rti, resp = respi, truth = truthi, err_tol = err_tol,
                     # limits:   vu,   vl,   a,      t0, w,  sv
                     lower = c(-Inf, -Inf, .01,       0, 0,   0),
                     upper = c( Inf,  Inf, Inf, min_rti, 1, Inf))
      res[["Convergence"]][(i-1)*ni+10*ninit_vals+j] <- temp[["convergence"]]
      res[["Objective"]][(i-1)*ni+10*ninit_vals+j] <- temp[["objective"]]
      res[["Iterations"]][(i-1)*ni+10*ninit_vals+j] <- temp[["iterations"]]
      res[["FuncEvals"]][(i-1)*ni+10*ninit_vals+j] <- temp[["evaluations"]][[1]]

      temp <- nlminb(init_vals[j, ], ll_fb_SWSE_14_5,
                     rt = rti, resp = respi, truth = truthi, err_tol = err_tol,
                     # limits:   vu,   vl,   a,      t0, w,  sv
                     lower = c(-Inf, -Inf, .01,       0, 0,   0),
                     upper = c( Inf,  Inf, Inf, min_rti, 1, Inf))
      res[["Convergence"]][(i-1)*ni+11*ninit_vals+j] <- temp[["convergence"]]
      res[["Objective"]][(i-1)*ni+11*ninit_vals+j] <- temp[["objective"]]
      res[["Iterations"]][(i-1)*ni+11*ninit_vals+j] <- temp[["iterations"]]
      res[["FuncEvals"]][(i-1)*ni+11*ninit_vals+j] <- temp[["evaluations"]][[1]]

      temp <- nlminb(init_vals[j, ], ll_fb_SWSE_17_6,
                     rt = rti, resp = respi, truth = truthi, err_tol = err_tol,
                     # limits:   vu,   vl,   a,      t0, w,  sv
                     lower = c(-Inf, -Inf, .01,       0, 0,   0),
                     upper = c( Inf,  Inf, Inf, min_rti, 1, Inf))
      res[["Convergence"]][(i-1)*ni+12*ninit_vals+j] <- temp[["convergence"]]
      res[["Objective"]][(i-1)*ni+12*ninit_vals+j] <- temp[["objective"]]
      res[["Iterations"]][(i-1)*ni+12*ninit_vals+j] <- temp[["iterations"]]
      res[["FuncEvals"]][(i-1)*ni+12*ninit_vals+j] <- temp[["evaluations"]][[1]]

      temp <- nlminb(init_vals[j, ], ll_fb_SWSE_14_6,
                     rt = rti, resp = respi, truth = truthi, err_tol = err_tol,
                     # limits:   vu,   vl,   a,      t0, w,  sv
                     lower = c(-Inf, -Inf, .01,       0, 0,   0),
                     upper = c( Inf,  Inf, Inf, min_rti, 1, Inf))
      res[["Convergence"]][(i-1)*ni+13*ninit_vals+j] <- temp[["convergence"]]
      res[["Objective"]][(i-1)*ni+13*ninit_vals+j] <- temp[["objective"]]
      res[["Iterations"]][(i-1)*ni+13*ninit_vals+j] <- temp[["iterations"]]
      res[["FuncEvals"]][(i-1)*ni+13*ninit_vals+j] <- temp[["evaluations"]][[1]]

      temp <- nlminb(init_vals[j, ], ll_fb_SWSE_17_7,
                     rt = rti, resp = respi, truth = truthi, err_tol = err_tol,
                     # limits:   vu,   vl,   a,      t0, w,  sv
                     lower = c(-Inf, -Inf, .01,       0, 0,   0),
                     upper = c( Inf,  Inf, Inf, min_rti, 1, Inf))
      res[["Convergence"]][(i-1)*ni+14*ninit_vals+j] <- temp[["convergence"]]
      res[["Objective"]][(i-1)*ni+14*ninit_vals+j] <- temp[["objective"]]
      res[["Iterations"]][(i-1)*ni+14*ninit_vals+j] <- temp[["iterations"]]
      res[["FuncEvals"]][(i-1)*ni+14*ninit_vals+j] <- temp[["evaluations"]][[1]]

      temp <- nlminb(init_vals[j, ], ll_fb_SWSE_14_7,
                     rt = rti, resp = respi, truth = truthi, err_tol = err_tol,
                     # limits:   vu,   vl,   a,      t0, w,  sv
                     lower = c(-Inf, -Inf, .01,       0, 0,   0),
                     upper = c( Inf,  Inf, Inf, min_rti, 1, Inf))
      res[["Convergence"]][(i-1)*ni+15*ninit_vals+j] <- temp[["convergence"]]
      res[["Objective"]][(i-1)*ni+15*ninit_vals+j] <- temp[["objective"]]
      res[["Iterations"]][(i-1)*ni+15*ninit_vals+j] <- temp[["iterations"]]
      res[["FuncEvals"]][(i-1)*ni+15*ninit_vals+j] <- temp[["evaluations"]][[1]]

      # microbenchmark
      mbm <- microbenchmark(
        fb_SWSE_17_0 = nlminb(init_vals[j,], ll_fb_SWSE_17_0, err_tol = err_tol,
                              rt = rti, resp = respi, truth = truthi,
                              # limits:   vu,   vl,   a,      t0, w,  sv
                              lower = c(-Inf, -Inf, .01,       0, 0,   0),
                              upper = c( Inf,  Inf, Inf, min_rti, 1, Inf)),
        fb_SWSE_14_0 = nlminb(init_vals[j,], ll_fb_SWSE_14_0, err_tol = err_tol,
                              rt = rti, resp = respi, truth = truthi,
                              # limits:   vu,   vl,   a,      t0, w,  sv
                              lower = c(-Inf, -Inf, .01,       0, 0,   0),
                              upper = c( Inf,  Inf, Inf, min_rti, 1, Inf)),
        fb_SWSE_17_1 = nlminb(init_vals[j,], ll_fb_SWSE_17_1, err_tol = err_tol,
                              rt = rti, resp = respi, truth = truthi,
                              # limits:   vu,   vl,   a,      t0, w,  sv
                              lower = c(-Inf, -Inf, .01,       0, 0,   0),
                              upper = c( Inf,  Inf, Inf, min_rti, 1, Inf)),
        fb_SWSE_14_1 = nlminb(init_vals[j,], ll_fb_SWSE_14_1, err_tol = err_tol,
                              rt = rti, resp = respi, truth = truthi,
                              # limits:   vu,   vl,   a,      t0, w,  sv
                              lower = c(-Inf, -Inf, .01,       0, 0,   0),
                              upper = c( Inf,  Inf, Inf, min_rti, 1, Inf)),
        fb_SWSE_17_2 = nlminb(init_vals[j,], ll_fb_SWSE_17_2, err_tol = err_tol,
                              rt = rti, resp = respi, truth = truthi,
                              # limits:   vu,   vl,   a,      t0, w,  sv
                              lower = c(-Inf, -Inf, .01,       0, 0,   0),
                              upper = c( Inf,  Inf, Inf, min_rti, 1, Inf)),
        fb_SWSE_14_2 = nlminb(init_vals[j,], ll_fb_SWSE_14_2, err_tol = err_tol,
                              rt = rti, resp = respi, truth = truthi,
                              # limits:   vu,   vl,   a,      t0, w,  sv
                              lower = c(-Inf, -Inf, .01,       0, 0,   0),
                              upper = c( Inf,  Inf, Inf, min_rti, 1, Inf)),
        fb_SWSE_17_3 = nlminb(init_vals[j,], ll_fb_SWSE_17_3, err_tol = err_tol,
                              rt = rti, resp = respi, truth = truthi,
                              # limits:   vu,   vl,   a,      t0, w,  sv
                              lower = c(-Inf, -Inf, .01,       0, 0,   0),
                              upper = c( Inf,  Inf, Inf, min_rti, 1, Inf)),
        fb_SWSE_14_3 = nlminb(init_vals[j,], ll_fb_SWSE_14_3, err_tol = err_tol,
                              rt = rti, resp = respi, truth = truthi,
                              # limits:   vu,   vl,   a,      t0, w,  sv
                              lower = c(-Inf, -Inf, .01,       0, 0,   0),
                              upper = c( Inf,  Inf, Inf, min_rti, 1, Inf)),
        fb_SWSE_17_4 = nlminb(init_vals[j,], ll_fb_SWSE_17_4, err_tol = err_tol,
                              rt = rti, resp = respi, truth = truthi,
                              # limits:   vu,   vl,   a,      t0, w,  sv
                              lower = c(-Inf, -Inf, .01,       0, 0,   0),
                              upper = c( Inf,  Inf, Inf, min_rti, 1, Inf)),
        fb_SWSE_14_4 = nlminb(init_vals[j,], ll_fb_SWSE_14_4, err_tol = err_tol,
                              rt = rti, resp = respi, truth = truthi,
                              # limits:   vu,   vl,   a,      t0, w,  sv
                              lower = c(-Inf, -Inf, .01,       0, 0,   0),
                              upper = c( Inf,  Inf, Inf, min_rti, 1, Inf)),
        fb_SWSE_17_5 = nlminb(init_vals[j,], ll_fb_SWSE_17_5, err_tol = err_tol,
                              rt = rti, resp = respi, truth = truthi,
                              # limits:   vu,   vl,   a,      t0, w,  sv
                              lower = c(-Inf, -Inf, .01,       0, 0,   0),
                              upper = c( Inf,  Inf, Inf, min_rti, 1, Inf)),
        fb_SWSE_14_5 = nlminb(init_vals[j,], ll_fb_SWSE_14_5, err_tol = err_tol,
                              rt = rti, resp = respi, truth = truthi,
                              # limits:   vu,   vl,   a,      t0, w,  sv
                              lower = c(-Inf, -Inf, .01,       0, 0,   0),
                              upper = c( Inf,  Inf, Inf, min_rti, 1, Inf)),
        fb_SWSE_17_6 = nlminb(init_vals[j,], ll_fb_SWSE_17_6, err_tol = err_tol,
                              rt = rti, resp = respi, truth = truthi,
                              # limits:   vu,   vl,   a,      t0, w,  sv
                              lower = c(-Inf, -Inf, .01,       0, 0,   0),
                              upper = c( Inf,  Inf, Inf, min_rti, 1, Inf)),
        fb_SWSE_14_6 = nlminb(init_vals[j,], ll_fb_SWSE_14_6, err_tol = err_tol,
                              rt = rti, resp = respi, truth = truthi,
                              # limits:   vu,   vl,   a,      t0, w,  sv
                              lower = c(-Inf, -Inf, .01,       0, 0,   0),
                              upper = c( Inf,  Inf, Inf, min_rti, 1, Inf)),
        fb_SWSE_17_7 = nlminb(init_vals[j,], ll_fb_SWSE_17_7, err_tol = err_tol,
                              rt = rti, resp = respi, truth = truthi,
                              # limits:   vu,   vl,   a,      t0, w,  sv
                              lower = c(-Inf, -Inf, .01,       0, 0,   0),
                              upper = c( Inf,  Inf, Inf, min_rti, 1, Inf)),
        fb_SWSE_14_7 = nlminb(init_vals[j,], ll_fb_SWSE_14_7, err_tol = err_tol,
                              rt = rti, resp = respi, truth = truthi,
                              # limits:   vu,   vl,   a,      t0, w,  sv
                              lower = c(-Inf, -Inf, .01,       0, 0,   0),
                              upper = c( Inf,  Inf, Inf, min_rti, 1, Inf)),
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



##### Benchmarks ###############################################################
data(med_dec, package = "fddm")
med_dec <- med_dec[which(med_dec[["rt"]] >= 0), ]
fit <- rt_fit(med_dec, id_idx = c(2,1), rt_idx = 8, response_idx = 7,
              truth_idx = 5, response_upper = "blast", err_tol = 1e-6,
              times = 5, unit = "ns")

save(fit, compress = "xz", compression_level = 9,
     file = paste0(save_dir, "delta_fit.Rds"))


##### Plot Results #############################################################
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

obj_diff_label <- function(y, df, col_name, mult = 1.15) {
  upper_limit <- max(df[[as.character(col_name)]]) * mult
  return(
    data.frame(
      y = 0.95 * upper_limit,
      label = paste(sum(y > 0, na.rm = TRUE))
    )
  )
}


Names <- c("fb_SWSE_17_0", "fb_SWSE_14_0", "fb_SWSE_17_1", "fb_SWSE_14_1",
           "fb_SWSE_17_2", "fb_SWSE_14_2", "fb_SWSE_17_3", "fb_SWSE_14_3",
           "fb_SWSE_17_4", "fb_SWSE_14_4", "fb_SWSE_17_5", "fb_SWSE_14_5",
           "fb_SWSE_17_6", "fb_SWSE_14_6", "fb_SWSE_17_7", "fb_SWSE_14_7")
Color <- c("#818679", "#cdcfc9", "#92c639", "#d3e8b0",
           "#8bb34d", "#d0e0b8", "#88ac53", "#d0deba",
           "#87a659", "#cfdbbd", "#869f60", "#cfd9bf",
           "#859966", "#ced6c2", "#83936c", "#cdd4c4")
Outline <- c("#818679", "#818679", "#92c639", "#92c639",
             "#8bb34d", "#8bb34d", "#88ac53", "#88ac53",
             "#87a659", "#87a659", "#869f60", "#869f60",
             "#859966", "#859966", "#83936c", "#83936c")
Shape <- c(21, 25)
Sizes <- c(0, 3, 3)
Stroke <- c(0, 1, 1)
Fills <- c("#ffffff00", "#ffffff00", "#80808099")



# uncomment the following line if loading pre-run benchmark data,
# will load into variable named "fit"
# load(paste0(save_dir, "delta_fit.Rds"))
fit <- fit_prep(fit)


# Benchmark Time Results
fit_mbm <- melt(fit, id.vars = c("Algorithm", "Convergence", "Obj_diff"),
                measure.vars = "BmTime", value.name = "BmTime")[,-4]

# Figure 5
fig_5 <- ggplot(fit_mbm, aes(x = factor(Algorithm, levels = Names),
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
  stat_summary(aes(y = Obj_diff, color = factor(Algorithm, levels = Names)),
               fun.data = obj_diff_label,
               fun.args = list(fit, "BmTime", 1.15),
               geom = "label",
               hjust = 0.5,
               vjust = 0.9) +
  scale_x_discrete(labels = c(bquote("0, " ~ S[17]), bquote("0, " ~ S[14]),
                              bquote("1, " ~ S[17]), bquote("1, " ~ S[14]),
                              bquote("2, " ~ S[17]), bquote("2, " ~ S[14]),
                              bquote("3, " ~ S[17]), bquote("3, " ~ S[14]),
                              bquote("4, " ~ S[17]), bquote("4, " ~ S[14]),
                              bquote("5, " ~ S[17]), bquote("5, " ~ S[14]),
                              bquote("6, " ~ S[17]), bquote("6, " ~ S[14]),
                              bquote("7, " ~ S[17]), bquote("7, " ~ S[14]))) +
  scale_color_manual(values = Outline, guide = FALSE) +
  scale_fill_manual(values = Color, guide = FALSE) +
  scale_shape_manual(values = Shape,
                     name = "Convergence Code",
                     breaks = c(0, 1),
                     labels = c("Success", "Failure")) +
  scale_size_manual(values = Sizes, guide = FALSE) +
  scale_discrete_manual(aesthetics = "stroke", values = Stroke, guide = FALSE) +
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
  labs(x = bquote(delta ~ ", and Summation Style"), y = "Time (microseconds)") +
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
ggsave(paste0(img_dir, "delta_fit_bm.png"),
       plot = fig_5, width = 16, height = 9)
