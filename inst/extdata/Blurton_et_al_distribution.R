# Distribution at lower barrier (Eq. 3 of the article)
# t: time (vector)
# nu: average drift
# eta2: variance of the drift distribution
# sigma2: variance of Wiener process
# a: upper barrier
# w: relative position of X(0) = z, w = z/a
# eps: required precision
#
G_0 = function(t=1.2, nu=0.1, eta2=0.01, sigma2=0.01, a=0.08, w=.375,
  eps=sqrt(.Machine$double.eps))
{
  nu = nu / sqrt(sigma2)
  a = a / sqrt(sigma2)
  eta2 = eta2 / sigma2
  sqt = sqrt(t)
  sqet = sqt * sqrt(1 + eta2*t)
  G = numeric(length(t))
  j = 0
  # Kendal added this to make eps (error tolerance) accurate inside the sum
  eps = eps / exp((-nu*nu*t - 2*nu*a*w + eta2*a*a*w*w)/2/(1 + eta2*t))
  repeat
  {
    rj = j*a + a*w
    logphi = dnorm(rj/sqt, log=TRUE)
    logM1 = logMill((rj - nu*t + eta2*(rj + a*w)*t) / sqet)
    logM2 = logMill((rj + nu*t + eta2*(rj - a*w)*t) / sqet)
    gj = exp(logphi + logM1) + exp(logphi + logM2)
    G = G + gj
    # Kendal is not sure why eps is only checked every 2 iterations
    if(all(gj < eps))
      return(exp((-nu*nu*t - 2*nu*a*w + eta2*a*a*w*w)/2/(1 + eta2*t)) * G)
    j = j + 1
    rj = j*a + a*(1-w)
    logphi = dnorm(rj/sqt, log=TRUE)
    logM1 = logMill((rj - nu*t + eta2*t*(rj + a*w)) / sqet)
    logM2 = logMill((rj + nu*t + eta2*t*(rj - a*w)) / sqet)
    gj = exp(logphi + logM1) + exp(logphi + logM2)
    G = G - gj
    j = j + 1
  }
}

# Distribution at upper barrier
#
G_a = function(t=1.2, nu=0.1, eta2=0.01, sigma2=0.01, a=0.08, w=.375,
  eps=sqrt(.Machine$double.eps))
{
  G_0(t, -nu, eta2, sigma2, a, 1-w, eps)
}

# log of Mill's ratio for the normal distribution
#
logMill = function(x) # log of Mill's ratio
{
  m = numeric(length(x))
  m[x >= 10000] = -log(x[x >= 10000]) # limiting case for x -> Inf
  m[x < 10000] = pnorm(x[x < 10000], lower=FALSE, log=TRUE) -
  dnorm(x[x < 10000], log=TRUE)
  m
}
