ezddm <- function(propCorrect, rtCorrectVariance_seconds, rtCorrectMean_seconds,
                  s = 0.1, nTrials) {

  # s is scaling parameter (defaults to 0.1 in Ratcliff's models)
  s2 <- s^2 # variance

  v <- as.numeric(NA)
  a <- as.numeric(NA)
  Ter <- as.numeric(NA)

  # Modify propCorrect (pC) because 0, 0.5, and 1 are all problematic values and
  # require edge correction. Original code added 0.00001 when pC = 0 or 0.5 and
  # changed pC = 1 - (1/2*nTrials) when pC = 1. The following transformation
  # adds eps^2 to pC if pC = 0, slightly reduces pC if pC = 0.5 or 1.
  # Note that in our use case:
  # 1) pC is always between 0 and 1, and nTrials is always given;
  # 2) we don't care about the absolute precision of the ezddm result as we are
  #    just looking to get decent initial values for the ddm fitting
  propCorrect <- (2 * nTrials * propCorrect) / (2 * nTrials + 1) +
                 1 / (4 * nTrials*nTrials)

  L <- stats::qlogis(propCorrect) # calculates logit
  x <- L * (L * propCorrect^2 - L * propCorrect + propCorrect - 0.5) /
       rtCorrectVariance_seconds
  v <- sign(propCorrect - 0.5) * s * x^(1/4) # drift rate
  a <- s2 * stats::qlogis(propCorrect) / v # threshold
  y <- -v * a / s2
  MDT <- (a / (2*v)) * (1 - exp(y)) / (1 + exp(y))
  Ter <- rtCorrectMean_seconds - MDT # non-decision time

  return(data.frame(a, v, Ter))
}
