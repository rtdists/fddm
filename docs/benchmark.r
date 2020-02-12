# R functions for running Microbenchmark tests on the response time models

# Utility packages and related local code
library("Rcpp")
library("microbenchmark")
sourceCpp("nearly_equal.cpp")
source("loglikelihood.r")

# Distribution packages and local code
library("rtdists")
library("RWiener")
source("BGK2014.r")





####################### Constant Drift Rate
rt_bm_cnst <- function(RT, V, A, W, t0=0.0001, eps=0.0, times=1000, unit="us"){
  nf <- 14 # number of functions being benchmarked
  nRT <- length(RT) # number of response times
  nV <- length(V) # number of drift rates
  nA <- length(A) # number of boundary separations
  nW <- length(W) # number of starting points
  fnames <- c("fs_eps_14", "fs_eps_17", "fs_BGK_14",
                 "fs_BGK_17", "fs_Nav_14", "fs_Nav_17",
                 "fb_BGK_Nav_14", "fb_BGK_Nav_17",
                 "fb_Nav_Nav_14", "fb_Nav_Nav_17", "fl_Nav",
                 "RWiener_R", "BGK2014_R", "RTDists_R")

  # Initialize the dataframe to contain the microbenchmark results
  mbm_res <- data.frame(matrix(ncol=times+5, nrow=nf*nRT*nV*nA*nW))
  colnames(mbm_res) <- c('RT', 'V', 'A', 'W', 'FuncName',
                         paste0("bm", seq_len(times)))
  start <- 1
  stop <- nf
  sl_time <- seq_len(times)+5

  # Loop through each combination of parameters and record microbenchmark results
  for (rt in 1:nRT) {
    for (v in 1:nV) {
      for (a in 1:nA) {
        for (w in 1:nW) {
          mbm <- microbenchmark(
            fs_eps_14 = fs_eps_2014(rt=RT[rt], v=V[v], a=A[a], w=W[w], t0=t0, eps=eps),
            fs_eps_17 = fs_eps_2017(rt=RT[rt], v=V[v], a=A[a], w=W[w], t0=t0, eps=eps),
            fs_BGK_14 = fs_BGK_2014(rt=RT[rt], v=V[v], a=A[a], w=W[w], t0=t0, eps=eps),
            fs_BGK_17 = fs_BGK_2017(rt=RT[rt], v=V[v], a=A[a], w=W[w], t0=t0, eps=eps),
            fs_Nav_14 = fs_Nav_2014(rt=RT[rt], v=V[v], a=A[a], w=W[w], t0=t0, eps=eps),
            fs_Nav_17 = fs_Nav_2017(rt=RT[rt], v=V[v], a=A[a], w=W[w], t0=t0, eps=eps),
            fb_BGK_Nav_14 = fb_BGK_Nav_14(rt=RT[rt], v=V[v], a=A[a], w=W[w], t0=t0, eps=eps),
            fb_BGK_Nav_17 = fb_BGK_Nav_17(rt=RT[rt], v=V[v], a=A[a], w=W[w], t0=t0, eps=eps),
            fb_Nav_Nav_14 = fb_Nav_Nav_14(rt=RT[rt], v=V[v], a=A[a], w=W[w], t0=t0, eps=eps),
            fb_Nav_Nav_17 = fb_Nav_Nav_17(rt=RT[rt], v=V[v], a=A[a], w=W[w], t0=t0, eps=eps),
            fl_Nav = fl_Nav(rt=RT[rt], v=V[v], a=A[a], w=W[w], t0=t0, eps=eps),
            RWiener_R = dwiener(RT[rt], resp="lower", delta=V[v], alpha=A[a], beta=W[w], tau=t0),
            BGK2014_R = fs14_R(t=RT[rt]-t0, v=V[v], a=A[a], w=W[w], eps=eps),
            RTDists_R = ddiffusion(RT[rt], "lower", v=V[v], a=A[a], z=W[w]*A[a], t0=t0),
            times = times, unit = unit)
          # add the rt, v, a, w values and function names to the dataframe
          mbm_res[start:stop, 1] <- rep(RT[rt], nf)
          mbm_res[start:stop, 2] <- rep(V[v]  , nf)
          mbm_res[start:stop, 3] <- rep(A[a]  , nf)
          mbm_res[start:stop, 4] <- rep(W[w]  , nf)
          mbm_res[start:stop, 5] <- fnames
          # add the microbenchmark results to the dataframe
          for (i in 1:nf) {
            bm <- subset(mbm, expr==fnames[i])
            mbm_res[start+i-1, sl_time] <- bm$time
          }
          # iterate start and stop values
          start = start + nf
          stop = stop + nf
        }
      }
    }
  }
  return(mbm_res)
}


rt_bm_cnst_vec <- function(RT, V, A, W, t0=0.0001, eps=0.0, times=1000, unit="us"){
  nf <- 14 # number of functions being benchmarked
  nV <- length(V) # number of drift rates
  nA <- length(A) # number of boundary separations
  nW <- length(W) # number of starting points
  resp <- rep("lower", length(RT)) # for RWiener
  fnames <- c("fs_eps_14", "fs_eps_17", "fs_BGK_14",
                 "fs_BGK_17", "fs_Nav_14", "fs_Nav_17",
                 "fb_BGK_Nav_14", "fb_BGK_Nav_17",
                 "fb_Nav_Nav_14", "fb_Nav_Nav_17", "fl_Nav",
                 "RWiener_R", "BGK2014_R", "RTDists_R")

  # Initialize the dataframe to contain the microbenchmark results
  mbm_res <- data.frame(matrix(ncol=times+4, nrow=nf*nV*nA*nW))
  colnames(mbm_res) <- c('V', 'A', 'W', 'FuncName',
                         paste0("bm", seq_len(times)))
  start <- 1
  stop <- nf
  sl_time <- seq_len(times)+4

  # Loop through each combination of parameters and record microbenchmark results
  for (v in 1:nV) {
    for (a in 1:nA) {
      for (w in 1:nW) {
        mbm <- microbenchmark(
          fs_eps_14 = fs_eps_2014(rt=RT, v=V[v], a=A[a], w=W[w], t0=t0, eps=eps),
          fs_eps_17 = fs_eps_2017(rt=RT, v=V[v], a=A[a], w=W[w], t0=t0, eps=eps),
          fs_BGK_14 = fs_BGK_2014(rt=RT, v=V[v], a=A[a], w=W[w], t0=t0, eps=eps),
          fs_BGK_17 = fs_BGK_2017(rt=RT, v=V[v], a=A[a], w=W[w], t0=t0, eps=eps),
          fs_Nav_14 = fs_Nav_2014(rt=RT, v=V[v], a=A[a], w=W[w], t0=t0, eps=eps),
          fs_Nav_17 = fs_Nav_2017(rt=RT, v=V[v], a=A[a], w=W[w], t0=t0, eps=eps),
          fb_BGK_Nav_14 = fb_BGK_Nav_14(rt=RT, v=V[v], a=A[a], w=W[w], t0=t0, eps=eps),
          fb_BGK_Nav_17 = fb_BGK_Nav_17(rt=RT, v=V[v], a=A[a], w=W[w], t0=t0, eps=eps),
          fb_Nav_Nav_14 = fb_Nav_Nav_14(rt=RT, v=V[v], a=A[a], w=W[w], t0=t0, eps=eps),
          fb_Nav_Nav_17 = fb_Nav_Nav_17(rt=RT, v=V[v], a=A[a], w=W[w], t0=t0, eps=eps),
          fl_Nav = fl_Nav(rt=RT, v=V[v], a=A[a], w=W[w], t0=t0, eps=eps),
          RWiener_R = dwiener(RT, resp=resp, delta=V[v], alpha=A[a], beta=W[w], tau=t0),
          BGK2014_R = fs14_R(t=RT-t0, v=V[v], a=A[a], w=W[w], eps=eps),
          RTDists_R = ddiffusion(RT, "lower", v=V[v], a=A[a], z=W[w]*A[a], t0=t0),
          times = times, unit = unit)
        # add the v, a, w values and function names to the dataframe
        mbm_res[start:stop, 1] <- rep(V[v]  , nf)
        mbm_res[start:stop, 2] <- rep(A[a]  , nf)
        mbm_res[start:stop, 3] <- rep(W[w]  , nf)
        mbm_res[start:stop, 4] <- fnames
        # add the microbenchmark results to the dataframe
        for (i in 1:nf) {
          bm <- subset(mbm, expr==fnames[i])
          mbm_res[start+i-1, sl_time] <- bm$time
        }
        # iterate start and stop values
        start = start + nf
        stop = stop + nf
      }
    }
  }
  return(mbm_res)
}





####################### Variable Drift Rate
rt_bm_vary <- function(RT, V, A, W, SV, t0=0.0001, eps=0.0, times=1000, unit="us") {
  nf <- 2 # number of functions being benchmarked
  nRT <- length(RT) # number of response times
  nV <- length(V) # number of drift rates
  nA <- length(A) # number of boundary separations
  nW <- length(W) # number of starting points
  nSV <- length(SV) # number of drift rate variances
  fnames <- c('RTDists_R', 'Vary_sv')

  # Initialize the dataframe to contain the microbenchmark results
  mbm_res <- data.frame(matrix(ncol=times+6, nrow=nf*nRT*nV*nA*nW*nSV))
  colnames(mbm_res) <- c('RT', 'V', 'A', 'W', 'SV', 'FuncName',
                         paste0("bm", as.character(seq_len(times))))
  start <- 1
  stop <- nf
  sl_time <- seq_len(times)+6

  # Loop through each combination of parameters and record microbenchmark results
  for (rt in 1:nRT) {
    for (v in 1:nV) {
      for (a in 1:nA) {
        for (w in 1:nW) {
          for (sv in 1:nSV) {
            mbm <- microbenchmark(
              RTDists_R = ddiffusion(RT[rt], "lower", v=V[v], a=A[a], z=W[w]*A[a],
                                     sv=SV[sv], t0=0),
              fs_vary = fs_vary(rt=RT[rt], v=V[v], a=A[a], w=W[w],
                                sv=SV[sv], eps=eps),
              times = times, unit = unit)
            # add the rt, v, a, w values and function names to the dataframe
            mbm_res[start:stop, 1] <- rep(RT[rt], nf)
            mbm_res[start:stop, 2] <- rep(V[v]  , nf)
            mbm_res[start:stop, 3] <- rep(A[a]  , nf)
            mbm_res[start:stop, 4] <- rep(W[w]  , nf)
            mbm_res[start:stop, 5] <- rep(SV[sv]  , nf)
            mbm_res[start:stop, 6] <- fnames
            # add the microbenchmark results to the dataframe
            for (i in 1:nf) {
              bm <- subset(mbm, expr==fnames[i])
              mbm_res[start+i-1, sl_time] <- bm$time
            }
            # iterate start and stop values
            start = start + nf
            stop = stop + nf
          }
        }
      }
    }
  }
  return(mbm_res)
}





####################### Wrapper for Benchmarks #################################
rt_benchmark <- function(RT, V, A, W, SV=0, t0=0.0001, eps=0.0,
                         RTvec=FALSE, times=1000, unit="us") {
  if (all(SV == 0)) { # constant drift rate
    if (RTvec) { # benchmark all rt's as a vector
      return(rt_bm_cnst_vec(RT=RT, V=V, A=A, W=W, t0=t0, eps=eps,
             times=times, unit=unit))
    } else { # benchmark each rt individually
      return(rt_bm_cnst(RT=RT, V=V, A=A, W=W, t0=t0, eps=eps,
             times=times, unit=unit))
    }
  } else { # variable drift rate
    return(rt_bm_vary(RT=RT, V=V, A=A, W=W, SV=SV, t0=t0, eps=eps,
           times=times, unit=unit))
  }
}




####################### Fitting to Real Data ###################################
rt_fit <- function(df, stvals=NULL, t0=0.0001,
                   eps=sqrt(.Machine$double.eps)) {

  # Preliminaries
  df <- subset(df, rt > t0) # remove any respnose times below t0

  inds <- sort(unique(df$ind))
  ninds <- max(length(inds), 1) # if inds is null, there is only one individual

  if (is.null(stvals)) {
    stvals <- list(#      v1,   v0,   a,   w
                   list(   0,    0,   1, 0.5), # defaults
                   list( 0.4, -0.4,   1, 0.5), # vary v1 & v0
                   list(-0.4,  0.4,   1, 0.5), # vary v1 & v0
                   list(   0,    0, 0.5, 0.5), # vary a
                   list(   0,    0, 2.5, 0.5), # vary a
                   list(   0,    0,   1, 0.2), # vary w
                   list(   0,    0,   1, 0.8)  # vary w
              )
  }
  nstvals <- length(stvals)
  stvals_mat <- matrix(unlist(stvals), ncol=4, nrow=7, byrow=TRUE)

  algo_names <- c(rep("fs_eps_17", nstvals), rep("fs_BGK_17", nstvals),
                  rep("fs_Nav_17", nstvals), rep("fl_Nav",    nstvals),
                  rep("RWiener_R", nstvals), rep("RTDists_R", nstvals))
  nalgos <- 6

  ni <- nalgos*nstvals

  # Initilize the result dataframe
  res <- data.frame(matrix(ncol = 11, nrow = ninds*nstvals*nalgos))
  colnames(res) <- c("ind", "Algorithm", "AlgoCalls", "initv1", "initv0",
                     "inita", "initw", "fitv1", "fitv0", "fita", "fitw")

  # Fill in what we can
  res[,2] <- algo_names # label algorithms
  res[,4] <- stvals_mat[,1] # label initial v1
  res[,5] <- stvals_mat[,2] # label initial v0
  res[,6] <- stvals_mat[,3] # label initial a
  res[,7] <- stvals_mat[,4] # label initial w

  # Loop through each individual and starting values
  for (i in 1:ninds) {
    dfi <- subset(df, ind == inds[i])
    res[((i-1)*ni+1):(i*ni), 1] <- inds[i] # label individuals
    rti <- dfi$rt
    respi <- dfi$response
    rrespi <- dfi$rresponse
    truthi <- dfi$truth
    for (j in 1:nstvals) {
      browser()
      temp <- nlminb(stvals[[j]], ll_fs_eps_17,
                     rt = rti, resp = respi, truth = truthi,
                     t0 = t0, eps = eps,
                     lower = c(-Inf, -Inf,   0, 0),
                     upper = c( Inf,  Inf, Inf, 1))
      res[(i-1)*ni+j,3] <- temp$evaluations[[1]]
      res[(i-1)*ni+j,8:11] <- temp$par

      temp <- nlminb(stvals[[j]], ll_fs_BGK_17,
                     rt = rti, resp = respi, truth = truthi,
                     t0 = t0, eps = eps,
                     lower = c(-Inf, -Inf,   0, 0),
                     upper = c( Inf,  Inf, Inf, 1))
      res[(i-1)*ni+nstvals+j,3] <- temp$evaluations[[1]]
      res[(i-1)*ni+nstvals+j,8:11] <- temp$par

      temp <- nlminb(stvals[[j]], ll_fs_Nav_17,
                     rt = rti, resp = respi, truth = truthi,
                     t0 = t0, eps = eps,
                     lower = c(-Inf, -Inf,   0, 0),
                     upper = c( Inf,  Inf, Inf, 1))
      res[(i-1)*ni+2*nstvals+j,3] <- temp$evaluations[[1]]
      res[(i-1)*ni+2*nstvals+j,8:11] <- temp$par

      temp <- nlminb(stvals[[j]], ll_fl_Nav,
                     rt = rti, resp = respi, truth = truthi,
                     t0 = t0, eps = eps,
                     lower = c(-Inf, -Inf,   0, 0),
                     upper = c( Inf,  Inf, Inf, 1))
      res[(i-1)*ni+3*nstvals+j,3] <- temp$evaluations[[1]]
      res[(i-1)*ni+3*nstvals+j,8:11] <- temp$par

      temp <- nlminb(stvals[[j]], ll_RWiener,
                     rt = rti, resp = rrespi, truth = truthi,
                     t0 = t0,
                     lower = c(-Inf, -Inf,   0, 0),
                     upper = c( Inf,  Inf, Inf, 1))
      res[(i-1)*ni+4*nstvals+j,3] <- temp$evaluations[[1]]
      res[(i-1)*ni+4*nstvals+j,8:11] <- temp$par

      temp <- nlminb(stvals[[j]], ll_RTDists,
                     rt = rti, resp = rrespi, truth = truthi,
                     t0 = t0,
                     lower = c(-Inf, -Inf,   0, 0),
                     upper = c( Inf,  Inf, Inf, 1))
      res[(i-1)*ni+5*nstvals+j,3] <- temp$evaluations[[1]]
      res[(i-1)*ni+5*nstvals+j,8:11] <- temp$par

    }
  }

  return(res)
}





####################### Accuracy Test for all Methods ##########################
rt_acc_test <- function(RT = c(0.1, seq(0.5, 3, by = 0.5), seq(4, 10, by = 1),
                               seq(12.5, 20, by = 2.5), seq(25, 30, b= 5)),
                        V = seq(-6, 6, by = 0.5),
                        A = seq(0.5, 5, by = 0.5),
                        W = seq(0.2, 0.8, by = 0.1),
                        t0 = 0.0001,
                        eps = sqrt(.Machine$double.eps),
                        Upper = FALSE){
  nf <- 14 # number of functions being benchmarked
  nRT <- length(RT) # number of response times
  nV <- length(V) # number of drift rates
  nA <- length(A) # number of boundary separations
  nW <- length(W) # number of starting points
  fnames <- c("fs_eps_14", "fs_eps_17", "fs_BGK_14",
              "fs_BGK_17", "fs_Nav_14", "fs_Nav_17",
              "fb_BGK_Nav_14", "fb_BGK_Nav_17",
              "fb_Nav_Nav_14", "fb_Nav_Nav_17", "fl_Nav",
              "RWiener_R", "BGK2014_R", "RTDists_R")
  resp0 = rep(0, nRT)
  rresp0 = rep("lower", nRT)

  if (!Upper) { # don't include upper density
    # Initialize the dataframe to contain the microbenchmark results
    res <- data.frame(matrix(ncol=8, nrow=nf*nRT*nV*nA*nW))
    colnames(res) <- c('rt', 'v', 'a', 'w', 'UpLo', 'FuncName', 'res', 'dif')
    start <- 1
    stop <- nf

    # Loop through each combination of parameters and record microbenchmark results
    for (rt in 1:nRT) {
      for (v in 1:nV) {
        for (a in 1:nA) {
          for (w in 1:nW) {
            # add the rt, v, a, w, upper/lower values, and function names to the dataframe
            res[start:stop, 1] <- rep(RT[rt], nf)
            res[start:stop, 2] <- rep(V[v]  , nf)
            res[start:stop, 3] <- rep(A[a]  , nf)
            res[start:stop, 4] <- rep(W[w]  , nf)
            res[start:stop, 5] <- rep("lower", nf)
            res[start:stop, 6] <- fnames

            # calculate "lower" density
            res[start,    7] <- fs_eps_2014(rt=RT[rt], resp=resp0[rt], v=V[v], a=A[a], w=W[w], t0=t0, eps=eps)
            res[start+1,  7] <- fs_eps_2017(rt=RT[rt], resp=resp0[rt], v=V[v], a=A[a], w=W[w], t0=t0, eps=eps)
            res[start+2,  7] <- fs_BGK_2014(rt=RT[rt], resp=resp0[rt], v=V[v], a=A[a], w=W[w], t0=t0, eps=eps)
            res[start+3,  7] <- fs_BGK_2017(rt=RT[rt], resp=resp0[rt], v=V[v], a=A[a], w=W[w], t0=t0, eps=eps)
            res[start+4,  7] <- fs_Nav_2014(rt=RT[rt], resp=resp0[rt], v=V[v], a=A[a], w=W[w], t0=t0, eps=eps)
            res[start+5,  7] <- fs_Nav_2017(rt=RT[rt], resp=resp0[rt], v=V[v], a=A[a], w=W[w], t0=t0, eps=eps)
            res[start+6,  7] <- fb_BGK_Nav_14(rt=RT[rt], resp=resp0[rt], v=V[v], a=A[a], w=W[w], t0=t0, eps=eps)
            res[start+7,  7] <- fb_BGK_Nav_17(rt=RT[rt], resp=resp0[rt], v=V[v], a=A[a], w=W[w], t0=t0, eps=eps)
            res[start+8,  7] <- fb_Nav_Nav_14(rt=RT[rt], resp=resp0[rt], v=V[v], a=A[a], w=W[w], t0=t0, eps=eps)
            res[start+9,  7] <- fb_Nav_Nav_17(rt=RT[rt], resp=resp0[rt], v=V[v], a=A[a], w=W[w], t0=t0, eps=eps)
            res[start+10, 7] <- fl_Nav(rt=RT[rt], resp=resp0[rt], v=V[v], a=A[a], w=W[w], t0=t0, eps=eps)
            res[start+11, 7] <- dwiener(RT[rt], resp=rresp0[rt], delta=V[v], alpha=A[a], beta=W[w], tau=t0)
            res[start+12, 7] <- fs14_R(t=RT[rt]-t0, v=V[v], a=A[a], w=W[w], eps=eps)
            res[start+13, 7] <- ddiffusion(RT[rt], rresp0[rt], v=V[v], a=A[a], z=W[w]*A[a], t0=t0)

            # determine accuracy of the results
            res[start:(start+13), 8] <- nearly_equal(res[start:(start+13), 7], eps)

            # iterate start and stop values
            start = start + nf
            stop = stop + nf
          }
        }
      }
    }

} else { # include upper density
    nf2 <- nf*2 # include upper/lower distributions
    resp1 = rep(1, nRT)
    rresp1 = rep("upper", nRT)

    # Initialize the dataframe to contain the microbenchmark results
    res <- data.frame(matrix(ncol=8, nrow=nf2*nRT*nV*nA*nW))
    colnames(res) <- c('rt', 'v', 'a', 'w', 'UpLo', 'FuncName', 'res', 'dif')
    start <- 1
    stop <- nf2

    # Loop through each combination of parameters and record microbenchmark results
    for (rt in 1:nRT) {
      for (v in 1:nV) {
        for (a in 1:nA) {
          for (w in 1:nW) {
            # add the rt, v, a, w, upper/lower values, and function names to the dataframe
            res[start:stop, 1] <- rep(RT[rt], nf2)
            res[start:stop, 2] <- rep(V[v]  , nf2)
            res[start:stop, 3] <- rep(A[a]  , nf2)
            res[start:stop, 4] <- rep(W[w]  , nf2)
            res[start:stop, 5] <- c(rep("lower", nf), rep("upper", nf))
            res[start:stop, 6] <- rep(fnames, 2)

            # calculate "lower" density
            res[start,    7] <- fs_eps_2014(rt=RT[rt], resp=resp0[rt], v=V[v], a=A[a], w=W[w], t0=t0, eps=eps)
            res[start+1,  7] <- fs_eps_2017(rt=RT[rt], resp=resp0[rt], v=V[v], a=A[a], w=W[w], t0=t0, eps=eps)
            res[start+2,  7] <- fs_BGK_2014(rt=RT[rt], resp=resp0[rt], v=V[v], a=A[a], w=W[w], t0=t0, eps=eps)
            res[start+3,  7] <- fs_BGK_2017(rt=RT[rt], resp=resp0[rt], v=V[v], a=A[a], w=W[w], t0=t0, eps=eps)
            res[start+4,  7] <- fs_Nav_2014(rt=RT[rt], resp=resp0[rt], v=V[v], a=A[a], w=W[w], t0=t0, eps=eps)
            res[start+5,  7] <- fs_Nav_2017(rt=RT[rt], resp=resp0[rt], v=V[v], a=A[a], w=W[w], t0=t0, eps=eps)
            res[start+6,  7] <- fb_BGK_Nav_14(rt=RT[rt], resp=resp0[rt], v=V[v], a=A[a], w=W[w], t0=t0, eps=eps)
            res[start+7,  7] <- fb_BGK_Nav_17(rt=RT[rt], resp=resp0[rt], v=V[v], a=A[a], w=W[w], t0=t0, eps=eps)
            res[start+8,  7] <- fb_Nav_Nav_14(rt=RT[rt], resp=resp0[rt], v=V[v], a=A[a], w=W[w], t0=t0, eps=eps)
            res[start+9,  7] <- fb_Nav_Nav_17(rt=RT[rt], resp=resp0[rt], v=V[v], a=A[a], w=W[w], t0=t0, eps=eps)
            res[start+10, 7] <- fl_Nav(rt=RT[rt], resp=resp0[rt], v=V[v], a=A[a], w=W[w], t0=t0, eps=eps)
            res[start+11, 7] <- dwiener(RT[rt], resp=rresp0[rt], delta=V[v], alpha=A[a], beta=W[w], tau=t0)
            res[start+12, 7] <- fs14_R(t=RT[rt]-t0, v=V[v], a=A[a], w=W[w], eps=eps)
            res[start+13, 7] <- ddiffusion(RT[rt], rresp0[rt], v=V[v], a=A[a], z=W[w]*A[a], t0=t0)

            # calculate "upper" density
            res[start+14, 7] <- fs_eps_2014(rt=RT[rt], resp=resp1[rt], v=V[v], a=A[a], w=W[w], t0=t0, eps=eps)
            res[start+15, 7] <- fs_eps_2017(rt=RT[rt], resp=resp1[rt], v=V[v], a=A[a], w=W[w], t0=t0, eps=eps)
            res[start+16, 7] <- fs_BGK_2014(rt=RT[rt], resp=resp1[rt], v=V[v], a=A[a], w=W[w], t0=t0, eps=eps)
            res[start+17, 7] <- fs_BGK_2017(rt=RT[rt], resp=resp1[rt], v=V[v], a=A[a], w=W[w], t0=t0, eps=eps)
            res[start+18, 7] <- fs_Nav_2014(rt=RT[rt], resp=resp1[rt], v=V[v], a=A[a], w=W[w], t0=t0, eps=eps)
            res[start+19, 7] <- fs_Nav_2017(rt=RT[rt], resp=resp1[rt], v=V[v], a=A[a], w=W[w], t0=t0, eps=eps)
            res[start+20, 7] <- fb_BGK_Nav_14(rt=RT[rt], resp=resp1[rt], v=V[v], a=A[a], w=W[w], t0=t0, eps=eps)
            res[start+21, 7] <- fb_BGK_Nav_17(rt=RT[rt], resp=resp1[rt], v=V[v], a=A[a], w=W[w], t0=t0, eps=eps)
            res[start+22, 7] <- fb_Nav_Nav_14(rt=RT[rt], resp=resp1[rt], v=V[v], a=A[a], w=W[w], t0=t0, eps=eps)
            res[start+23, 7] <- fb_Nav_Nav_17(rt=RT[rt], resp=resp1[rt], v=V[v], a=A[a], w=W[w], t0=t0, eps=eps)
            res[start+24, 7] <- fl_Nav(rt=RT[rt], resp=resp1[rt], v=V[v], a=A[a], w=W[w], t0=t0, eps=eps)
            res[start+25, 7] <- dwiener(RT[rt], resp=rresp1[rt], delta=V[v], alpha=A[a], beta=W[w], tau=t0)
            res[start+26, 7] <- fs14_R(t=RT[rt]-t0, v=-V[v], a=A[a], w=1-W[w], eps=eps) # use upper distribution parameters
            res[start+27, 7] <- ddiffusion(RT[rt], rresp1[rt], v=V[v], a=A[a], z=W[w]*A[a], t0=t0)

            # determine accuracy of the results
            res[start:(start+13), 8] <- nearly_equal(res[start:(start+13), 7], eps)
            res[(start+14):(start+27), 8] <- nearly_equal(res[(start+14):(start+27), 7], eps)

            # iterate start and stop values
            start = start + nf2
            stop = stop + nf2
          }
        }
      }
    }
  }

  return(res)
}






######### Visualizations of these results can be found in Code/Paper Images
