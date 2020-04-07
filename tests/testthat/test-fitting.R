context("Accuracy of fitting results across methods and established packages")
source("benchmark_testing/fitting_functions.r")


# Import data
data <- read.csv("inst/extdata/medical_dm.csv", stringsAsFactors = FALSE)

# Run/Import fitting
fit <- rt_fit(data, ind_idx = 1, rt_idx = 8, response_idx = 7, truth_idx = 5,
              response_upper = "blast")
fit <- readRDS("benchmark_testing/Results/fit.Rds")
eps <- 1e-4

# Calculate differences and means
fit$v1_dif <- rep(0, nrow(fit))
fit$v0_dif <- rep(0, nrow(fit))
fit$a_dif <- rep(0, nrow(fit))
fit$w_dif <- rep(0, nrow(fit))

inds <- unique(fit$ind)
ninds <- length(inds)
algos <- unique(fit$Algorithm)
nalgos <- length(algos)

avg_df <- data.frame(matrix(ncol = 6, nrow = ninds*nalgos))
colnames(avg_df) <- c('ind', 'Algorithm', 'v1_avg', 'v0_avg', 'a_avg', 'w_avg')
avg_df$Algorithm <- rep(algos, ninds)

start <- 1
stop <- nrow(subset(fit, ind == inds[1] & Algorithm == algos[1]))
for (i in 1:ninds) {
  avg_df$ind[((i-1)*nalgos+1):(i*nalgos)] <- rep(inds[i], nalgos) # label ind
  for (j in 1:nalgos) {
    tmp <- subset(fit, ind == inds[i] & Algorithm == algos[j])[,8:11]
    avg_df[(i-1)*nalgos+j,3:6] <- colMeans(tmp) # collect means
    for (k in 1:(stop-1)) {
      fit[start+k, 12:15] <- abs(fit[start, 8:11] - fit[(start+k), 8:11])
    }
    start = start + stop
  }
}


# Subset
Foster <- subset(fit, Algorithm == "fs_Fos_17")
Kesselmeier <- subset(fit, Algorithm == "fb_Kes_17")
Navarro <- subset(fit, Algorithm == "fb_Nav_17")
rtdists <- subset(fit, Algorithm == "rtdists")
RWiener <- subset(fit, Algorithm == "RWiener")
Kesselmeier_R <- subset(fit, Algorithm == "Kesselmeier")

Foster_avg <- subset(avg_df, Algorithm == "fs_Fos_17")
Kesselmeier_avg <- subset(avg_df, Algorithm == "fb_Kes_17")
Navarro_avg <- subset(avg_df, Algorithm == "fb_Nav_17")
rtdists_avg <- subset(avg_df, Algorithm == "rtdists")
RWiener_avg <- subset(avg_df, Algorithm == "RWiener")
Kesselmeier_R_avg <- subset(avg_df, Algorithm == "Kesselmeier")


# Test accuracy for fitted parameters: v_upper, v_lower, a, w
test_that("Accuracy within each method, using different starting values", {
  expect_true(all(Foster[,12:15] < eps))
  expect_true(all(Kesselmeier[,12:15] < eps))
  expect_true(all(Navarro[,12:15] < eps))
  expect_true(all(rtdists[,12:15] < eps))
  expect_true(all(RWiener[,12:15] < eps))
  # expect_true(all(Kesselmeier_R[,12:15] < eps))
})

test_that("Accuracy across different methods", {
  expect_true(all(abs(Foster_avg[,-(1:2)] - Kesselmeier_avg[,-(1:2)]) < eps))
  expect_true(all(abs(Foster_avg[,-(1:2)] - Navarro_avg[,-(1:2)]) < eps))
  expect_true(all(abs(Foster_avg[,-(1:2)] - rtdists_avg[,-(1:2)]) < eps))
  expect_true(all(abs(Foster_avg[,-(1:2)] - RWiener_avg[,-(1:2)]) < eps))
  # expect_true(all(abs(Foster_avg[,-(1:2)] - Kesselmeier_R_avg[,-(1:2)]) < eps))
})
