context("Comparing saved with new fits (validity vignette)")

test_that("Fits in validity vignette", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("rtdists")
  library("rtdists")

  ll_fb_SWSE_17 <- function(pars, rt, resp, truth, err_tol) {
    v <- numeric(length(rt))
    v[truth == "upper"] <- pars[[1]]
    v[truth == "lower"] <- pars[[2]]
    dens <- dfddm(rt = rt, response = resp, a = pars[[3]], v = v,
                  t0 = pars[[4]], w = pars[[5]], sv = pars[[6]], log = TRUE,
                  n_terms_small = "SWSE", summation_small = "2017",
                  scale = "both", err_tol = err_tol)
    return( ifelse(any(!is.finite(dens)), 1e6, -sum(dens)) )
  }

  ll_fb_Gon_17 <- function(pars, rt, resp, truth, err_tol) {
    v <- numeric(length(rt))
    v[truth == "upper"] <- pars[[1]]
    v[truth == "lower"] <- pars[[2]]
    dens <- dfddm(rt = rt, response = resp, a = pars[[3]], v = v,
                  t0 = pars[[4]], w = pars[[5]], sv = pars[[6]], log = TRUE,
                  n_terms_small = "Gondan", summation_small = "2017",
                  scale = "both", err_tol = err_tol)
    return( ifelse(any(!is.finite(dens)), 1e6, -sum(dens)) )
  }

  ll_fb_Nav_17 <- function(pars, rt, resp, truth, err_tol) {
    v <- numeric(length(rt))
    v[truth == "upper"] <- pars[[1]]
    v[truth == "lower"] <- pars[[2]]
    dens <- dfddm(rt = rt, response = resp, a = pars[[3]], v = v,
                  t0 = pars[[4]], w = pars[[5]], sv = pars[[6]], log = TRUE,
                  n_terms_small = "Navarro", summation_small = "2017",
                  scale = "both", err_tol = err_tol)
    return( ifelse(any(!is.finite(dens)), 1e6, -sum(dens)) )
  }

  ll_RTDists <- function(pars, rt, resp, truth) {
    v <- numeric(length(rt))
    v[truth == "upper"] <- pars[[1]]
    v[truth == "lower"] <- pars[[2]]
    dens <- log(ddiffusion(rt, resp, a = pars[[3]], v = v, t0 = pars[[4]],
                           z = pars[[5]]*pars[[3]], sv = pars[[6]]))
    return( ifelse(any(!is.finite(dens)), 1e6, -sum(dens)) )
  }

  rt_fit <- function(data, id_idx = NULL, rt_idx = NULL, response_idx = NULL,
                     truth_idx = NULL, response_upper = NULL, err_tol = 1e-6) {

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

    init_vals <- data.frame(vu = c( 0,  10, -.5,  0,  0,  0,  0,  0,  0,   0,  0),
                            vl = c( 0, -10,  .5,  0,  0,  0,  0,  0,  0,   0,  0),
                            a  = c( 1,   1,   1, .5,  5,  1,  1,  1,  1,   1,  1),
                            t0 = c( 0,   0,   0,  0,  0,  0,  0,  0,  0,   0,  0),
                            w  = c(.5,  .5,  .5, .5, .5, .5, .5, .2, .8,  .5, .5),
                            sv = c( 1,   1,   1,  1,  1,  1,  1,  1,  1, .05,  5))
    ninit_vals <- nrow(init_vals)

    algo_names <- c("fb_SWSE_17", "fb_Gon_17", "fb_Nav_17", "rtdists")
    nalgos <- length(algo_names)
    ni <- nalgos*ninit_vals

    # Initilize the result dataframe
    cnames <- c("ID", "Algorithm", "Convergence", "Objective",
                "vu_init", "vl_init", "a_init", "t0_init", "w_init", "sv_init",
                "vu_fit", "vl_fit", "a_fit", "t0_fit", "w_fit", "sv_fit")
    res <- data.frame(matrix(ncol = length(cnames), nrow = nids*ninit_vals*nalgos))
    colnames(res) <- cnames

    # label the result dataframe
    res[["ID"]] <- rep(ids, each = ni) # label individuals
    res[["Algorithm"]] <- rep(algo_names, each = ninit_vals) # label algorithms
    res[["vu_init"]] <- init_vals[["vu"]] # label initial vu
    res[["vl_init"]] <- init_vals[["vl"]] # label initial vl
    res[["a_init"]]  <- init_vals[["a"]]  # label initial a
    res[["w_init"]]  <- init_vals[["w"]]  # label initial w
    res[["sv_init"]] <- init_vals[["sv"]] # label initial sv

    # Loop through each individual and starting values
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

      # label the result dataframe
      res[["t0_init"]][((i-1)*ni+1):(i*ni)] <- init_vals[["t0"]] # label initial t0

      # loop through all of the starting values
      for (j in 1:ninit_vals) {
        temp <- nlminb(init_vals[j, ], ll_fb_SWSE_17,
                       rt = rti, resp = respi, truth = truthi, err_tol = err_tol,
                       # limits:   vu,   vl,   a,      t0, w,  sv
                       lower = c(-Inf, -Inf, .01,       0, 0,   0),
                       upper = c( Inf,  Inf, Inf, min_rti, 1, Inf))
        res[["Convergence"]][(i-1)*ni+0*ninit_vals+j] <- temp[["convergence"]]
        res[["Objective"]][(i-1)*ni+0*ninit_vals+j] <- temp[["objective"]]
        res[(i-1)*ni+0*ninit_vals+j, 11:16] <- temp[["par"]]

        temp <- nlminb(init_vals[j, ], ll_fb_Gon_17,
                       rt = rti, resp = respi, truth = truthi, err_tol = err_tol,
                       # limits:   vu,   vl,   a,      t0, w,  sv
                       lower = c(-Inf, -Inf, .01,       0, 0,   0),
                       upper = c( Inf,  Inf, Inf, min_rti, 1, Inf))
        res[["Convergence"]][(i-1)*ni+1*ninit_vals+j] <- temp[["convergence"]]
        res[["Objective"]][(i-1)*ni+1*ninit_vals+j] <- temp[["objective"]]
        res[(i-1)*ni+1*ninit_vals+j, 11:16] <- temp[["par"]]

        temp <- nlminb(init_vals[j, ], ll_fb_Nav_17,
                       rt = rti, resp = respi, truth = truthi, err_tol = err_tol,
                       # limits:   vu,   vl,   a,      t0, w,  sv
                       lower = c(-Inf, -Inf, .01,       0, 0,   0),
                       upper = c( Inf,  Inf, Inf, min_rti, 1, Inf))
        res[["Convergence"]][(i-1)*ni+2*ninit_vals+j] <- temp[["convergence"]]
        res[["Objective"]][(i-1)*ni+2*ninit_vals+j] <- temp[["objective"]]
        res[(i-1)*ni+2*ninit_vals+j, 11:16] <- temp[["par"]]

        temp <- nlminb(init_vals[j, ], ll_RTDists,
                       rt = rti, resp = respi, truth = truthi,
                       # limits:   vu,   vl,   a,      t0, w,  sv
                       lower = c(-Inf, -Inf, .01,       0, 0,   0),
                       upper = c( Inf,  Inf, Inf, min_rti, 1, Inf))
        res[["Convergence"]][(i-1)*ni+3*ninit_vals+j] <- temp[["convergence"]]
        res[["Objective"]][(i-1)*ni+3*ninit_vals+j] <- temp[["objective"]]
        res[(i-1)*ni+3*ninit_vals+j, 11:16] <- temp[["par"]]
      }
    }
    return(res)
  }
  data(med_dec, package = "fddm")
  med_dec <- med_dec[which(med_dec[["rt"]] >= 0),]
  newfit <- rt_fit(med_dec, id_idx = c(2,1), rt_idx = 8, response_idx = 7,
                   truth_idx = 5, response_upper = "blast", err_tol = 1e-6)
  load(system.file("extdata", "dfddm_density", "valid_fit.Rds", package = "fddm", mustWork = TRUE))
  expect_equal(newfit[["Objective"]], fit[["Objective"]], tolerance = 0.01)
})
