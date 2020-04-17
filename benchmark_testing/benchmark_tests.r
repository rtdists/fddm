library("devtools")
load_all(recompile = TRUE)
library("fddm")
source("benchmark_testing/benchmark_functions.r")


# Set parameter exploration vectors
RT_0_30 = c(0.001, 0.005, 0.01, 0.05, 0.1,
            seq(0.5, 3, by = 0.5),
            seq(4, 10, by = 1),
            seq(12.5, 20, by = 2.5),
            seq(25, 30, by = 5))
RT_0_3 = c(0.001, 0.005, 0.01, 0.05, 0.1,
           seq(0.5, 3, by = 0.5))
RT_4_10 = seq(4, 10, by = 1)
RT_12_30 = c(seq(12.5, 20, by = 2.5),
             seq(25, 30, by = 5))
A = c(0.25, seq(0.5, 5, by = 0.5))
V = seq(-6, 6, by = 0.5)
t0 = 1e-4 # must be nonzero for RWiener
W = seq(0.2, 0.8, by = 0.1)
err_tol = 1e-6 # this is the setting from rtdists


# Run benchmark tests
sim <- rt_benchmark(RT = c(RT_0_3, RT_4_10), resp = 0, V = V, A = A, W = W, t0 = t0,
                    err_tol = err_tol, rt_as_vec = FALSE,
                    times = 1000, unit = "us")
saveRDS(sim, file = "benchmark_testing/Results/sim_0-10_1000w_TEST3.Rds")
sim <- readRDS("benchmark_testing/Results/sim_0-3_10000.Rds")

sim_vec <- rt_benchmark(RT = RT_4_10, resp = 0, V = V, A = A, W = W, t0 = t0,
                        err_tol = err_tol, rt_as_vec = TRUE,
                        times = 10000, unit = "us")
saveRDS(sim_vec, file = "benchmark_testing/Results/vec_4-10_10000.Rds")
sim_vec <- readRDS("benchmark_testing/Results/vec_0-3_100000.Rds")