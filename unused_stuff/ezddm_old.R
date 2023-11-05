
ezddm <- function(propCorrect, rtCorrectVariance_seconds, rtCorrectMean_seconds, s = 0.1, nTrials = NULL) {

    # s is scaling parameter (defaults to 0.1 in Ratcliff's models)
    s2 <- s^2 # variance

    v <- as.numeric(NA)
    a <- as.numeric(NA)
    Ter <- as.numeric(NA)

    # if propCorrect equals 0, 0.5, or 1, this method will not work, and an edge correction is required
    if (propCorrect %in% c(0, 0.5, 1)) {

        if (propCorrect == 0) {
            propCorrect <- propCorrect + 0.00001
        } else if (propCorrect == 0.5) {
            #cat("Oops, propCorrect == 0.5 (chance performance; drift will be close to 0). Added 0.00001 to propCorrect.\n")
            propCorrect <- propCorrect + 0.00001
        } else if (propCorrect == 1) {
            if (!is.null(nTrials)) {
                #cat("Oops, propCorrect == 1. Applied edge correction.\n")
                propCorrect <- 1 - (1 / (2 * nTrials))
            } else {
                #cat("Oops, propCorrect == 1. Edge correction required. Provide number of trials (nTrials).\n")
            }
        }

    }

    if (propCorrect != 1) {

        L <- stats::qlogis(propCorrect) # calculates logit
        x <- L * (L * propCorrect^2 - L * propCorrect + propCorrect - 0.5) / rtCorrectVariance_seconds
        v <- sign(propCorrect - 0.5) * s * x^(1/4) # drift rate
        a <- s2 * stats::qlogis(propCorrect)/v # threshold
        y <- -v*a/s2
        MDT <- (a/(2*v)) * (1-exp(y))/(1 + exp(y))
        Ter <- rtCorrectMean_seconds - MDT # non-decision time

    }

    return(data.frame(a, v, Ter))
}





edgeCorrect <- function(n) {
    return(1 - (1 / (2 * n))) # n: number of observations
}
