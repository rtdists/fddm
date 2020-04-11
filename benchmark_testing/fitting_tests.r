library("devtools")
load_all(recompile = TRUE)
library("fddm")
source("benchmark_testing/fitting_functions.r")


# Import data
data <- read.csv("inst/extdata/medical_dm.csv", stringsAsFactors = FALSE)

# Run fitting
fit <- rt_fit(data, ind_idx = 1, rt_idx = 8, response_idx = 7, truth_idx = 5,
              response_upper = "blast")
saveRDS(fit, file = "benchmark_testing/Results/fit.Rds")
fit <- readRDS("benchmark_testing/Results/fit.Rds")
