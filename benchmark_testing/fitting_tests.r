library("devtools")
load_all(recompile = TRUE)
library("fddm")
source("benchmark_testing/fitting_functions.r")


# Import data
data(med_dec, package = "fddm")
data <- subset(med_dec, id == unique(med_dec$id)[1])
data <- subset(data, rt > t0)
# Run fitting
fit <- rt_fit(dat, ind_idx = 1, rt_idx = 8, response_idx = 7, truth_idx = 5,
              response_upper = "blast", times = 10, unit = "us")
saveRDS(fit, file = "benchmark_testing/Results/fit_10000.Rds")
fit <- readRDS("benchmark_testing/Results/fit_10000.Rds")




ll_fs_Fos_17 <- function(pars, rt, resp, truth, t0, err_tol) {
  rt1 <- rt[truth == 1]
  rt0 <- rt[truth == 0]
  resp1 <- resp[truth == 1]
  resp0 <- resp[truth == 0]

  # the truth is "upper" so use v1
  dens1 <- dfddm(rt = rt1, response = resp1, a = pars[3],
                 v = pars[1], t0 = t0, w = pars[4], sv = pars[5],
                 log = FALSE, n_terms_small = "Foster",
                 summation_small = "2017", scale = "small",
                 err_tol = err_tol)
  # the truth is "lower" so use v0
  dens0 <- dfddm(rt = rt0, response = resp0, a = pars[3],
                 v = pars[2], t0 = t0, w = pars[4], sv = pars[5],
                 log = FALSE, n_terms_small = "Foster",
                 summation_small = "2017", scale = "small",
                 err_tol = err_tol)

  densities <- c(dens1, dens0)
  if (any(!is.finite(densities)) | any(densities <= 0)) return(1e6)
  return(-sum(log(densities)))
}

any(c(seq_len(6), NA) == 0)
seq_len(length(dens1))[dens1 == 0]
dens1[124]
rt1[124]
resp1[124]

ind_idx = 1
rt_idx = 8
response_idx = 7
truth_idx = 5
response_upper = "blast"

nr <- nrow(data)
df <- data.frame(ind = integer(nr),
                 rt = double(nr),
                 response = integer(nr),
                 rresponse = character(nr),
                 truth = integer(nr))

df$ind <- as.integer(data[,ind_idx])
df$rt <- as.double(data[,rt_idx])

df$response <- as.integer(0)
df$response[data[,response_idx] == response_upper] <- as.integer(1)
df$rresponse <- "lower"
df$rresponse[data[,response_idx] == response_upper] <- "upper"

df$truth <- as.integer(0)
df$truth[data[,truth_idx] == response_upper] <- as.integer(1)
head(df)
pars <- c(   0,    0,   1, 0.5, 0.1)
rt <- df$rt
resp <- df$response
truth <- df$truth
t0 <- 1e-4
err_tol <- 1e-6



ll_fs_Fos_17(pars, rt, resp, truth, t0, err_tol)
nlminb(pars, ll_fs_Fos_17, rt = rt, resp = resp, truth = truth, t0 = t0, err_tol = err_tol,
  lower = c(-Inf, -Inf,   0, 0, 0.01),
  upper = c( Inf,  Inf, Inf, 1,  Inf))$par
