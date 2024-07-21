devtools::load_all(recompile = TRUE)
devtools::load_all()
library("fddm")

data <- read.csv("unused_stuff/ticket_price_data.csv")
data <- data[ !is.na(data[["response"]]), ]
data <- data[ data[["rt"]] >= 0.25 & data[["rt"]] <= 3, ]
data$is_first <- factor(ifelse(data$pos == 1, "first", "other"),
                        c("first", "other"))
data$response = data$response + 1
data$price_p <- pnorm(data$prices, 180, 20)
data$price_q <- (data$prices - 180) / 20

# data <- data[-c(1, 2, 3, 4, 5),]

m1 <- ddm(drift = rt + response ~ price_q,
          boundary = ~1, # ~pos:is_first,
          ndt = ~pos:is_first,
          bias = ~1, # ~pos:is_first,
          use_gradient = FALSE,
          optim = "optim",
          # optim = "nlminb",
          args_optim = list(
                            # init = c(-0.8319831, -0.7425308, 2.0484642,
                            #          -0.3685364, -0.1007685, 0.0500000,
                            #          0.0500000, 0.0500000, 0.5000000,
                            #          0.0000000, 0.0000000),
                            # control = list(eval.max = 5000,
                            #                iter.max = 5000)),
                            # method = "L-BFGS-B",
                            control = list(maxit = 5000)),
          data = data)

m1$optim_info$args_optim
m1$optim_info$value


### nlminb wo grad
# $par
# [1] -1.01780037 -0.80081309  2.16435476  0.13706442 -0.13706432  0.01129356
# [7]  0.50465685

# $objective
# [1] 90898.99

# $convergence
# [1] 0

# $iterations
# [1] 171

# $evaluations
# function gradient 
#      293     2086 

# $message
# [1] "relative convergence (4)"


### nlminb with grad
# $par
# [1] -0.607823962 -0.006275201  1.439077998  0.130570055 -0.004404149
# [6]  0.011942994  0.577434662

# $objective
# [1] 138347.7

# $convergence
# [1] 1

# $iterations
# [1] 36

# $evaluations
# function gradient 
#       81       36 

# $message
# [1] "false convergence (8)"


### optim wo grad
# $par
# [1] -1.414222399 -0.989105740  1.999470445  0.214294767  0.028366527
# [6]  0.001610306  0.552966568

# $value
# [1] 83609.88

# $counts
# function gradient 
#      876       NA 

# $convergence
# [1] 0

# $message
# NULL


### optim with grad
# $par
# [1] -1.788494389  0.001542208  1.336193716  0.147650614  0.002321144
# [6]  0.010228478  0.583601717

# $value
# [1] 197762.8

# $counts
# function gradient 
#       91       91 

# $convergence
# [1] 52

# $message
# [1] "ERROR: ABNORMAL_TERMINATION_IN_LNSRCH"




rt <- data[["rt"]]
resp <- data[["response"]]
v <- -0.0076043
a <- 1.12276
t0 <- 0.24994
w <- 0.497403
sv <- 0
# dens <- dfddm(rt = rt, response = resp, v = v, a = a, t0 = t0, w = w, sv = sv)
# rt0 <- rt[which(dens == 0)]
# resp0 <- resp[which(dens == 0)]
# dens0 <- dfddm(rt = rt0, response = resp0, v = v, a = a, t0 = t0, w = w, sv = sv)
dfddm(rt = 6e-5, resp = 1, v = v, a = a, t0 = 0, w = w, sv = sv, log = TRUE)
dfddm(rt = 10^(-2:-6), resp = 1, v = v, a = a, t0 = 0, w = w, sv = sv)





devtools::load_all()
devtools::load_all(recompile = TRUE)
library("fddm")


rt <- 1.256
resp <- 1
v <- -51.5559
# v <- -324.457
a <- 3592.06
t0 <- 0.227905
w <- 0.528234
sv <- 0
dfddm(rt = rt, resp = 1, v = v, a = a, t0 = t0, w = w, sv = sv)
t <- rt - t0
t / (a*a)
a * exp(-v * a * w - v * v * t / 2) / (t * sqrt(pi) * sqrt(t))
exp(10^(2:6))






sum(m1$compiled_model$sv <= 0)
sum(m1$compiled_model$sv < 0)
any(is.infinite(m1$compiled_model$sv))

m1$compiled_model$calculate_loglik(m1$compiled_model$coefficients)

fitted_coefs <- c(-0.83198308, -0.74253081, 2.04846421, -0.36853645,
                  -0.10076849, 0.13250274, 0.39634294, 0.03000149,
                  0.50000000, 0.00000000, 0.00000000)

# need to edit R/ddm_methods.R funcs for printing:
# Error in x$link$drift : $ operator is invalid for atomic vectors
