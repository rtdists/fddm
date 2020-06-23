#
# Online supplemental material for ``Even faster and more accurate
# first-passage time densities and distributions for the Wiener diffusion
# model'' by Gondan, Blurton, and Kesselmeier.
#

#
# log of Mill's ratio for the normal distribution
#
logMill = function(x) # log of Mills ratio
{
    m = numeric(length(x))
    m[x >= 10000] = -log(x[x >= 10000]) # limiting case for x -> Inf
    m[x < 10000]  = pnorm(x[x < 10000], lower=FALSE, log=TRUE) -
        dnorm(x[x < 10000], log=TRUE)
    return(m)
}

#
# Number of terms required for the distribution
#
Ks = function(t, v, a, w, eps)
{
    K1  = (abs(v)*t - a*w)/2/a
    arg = pmax(0, pmin(1, exp(v*a*w + v*v*t/2 + log(eps))/2))
    K2  = -sqrt(t)/2/a * qnorm(arg)
    return(ceiling(max(K1, K1 + K2)))
}

#
# Distribution at lower barrier - small time representation
#
# t: time (vector)
# v: drift
# a: upper barrier
# w: relative position of X(0) = z, w = z/a
# eps: required precision
#
Fs = function(t, v, a, w, eps=sqrt(.Machine$double.eps))
{
	K = Ks(t, v, a, w, eps)
	F = numeric(length(t))
	sqt = sqrt(t)
	for(k in K:0)
	{
		rj = 2*k*a + a*w
		dj = -v*a*w - v*v*t/2 + dnorm(rj/sqt, log=TRUE)
		pos1 = dj + logMill((rj-v*t)/sqt)
		pos2 = dj + logMill((rj+v*t)/sqt)
		rj = (2*k+1)*a + a*(1-w)
		dj = -v*a*w - v*v*t/2 + dnorm(rj/sqt, log=TRUE)
		neg1 = dj + logMill((rj-v*t)/sqt)
		neg2 = dj + logMill((rj+v*t)/sqt)
		F = exp(pos1) + exp(pos2) - exp(neg1) - exp(neg2) + F
	}
	return(F)
}

#
# Number of terms for the density
#
ks = function(t, w, eps)
{
	K1 = K2 = (sqrt(2*t) - w)/2
	u_eps = pmin(-1, log(2*pi*t*t*eps*eps)) # Safe bound so that
	arg = -t * (u_eps - sqrt(-2*u_eps - 2)) # sqrt(x) with x > 0
	K2[arg > 0] = 1/2 * sqrt(arg) - w/2
	return(ceiling(max(K1, K2)))
}

#
# Density at the lower barrier - one-parameter form
#
fsw = function(t, w, eps)
{
    K = ks(t, w, eps)
    f = numeric(length(t))
    if(K > 0) for(k in K:1)
        f = (w+2*k) * exp(-(w+2*k) * (w+2*k)/2/t) +
            (w-2*k) * exp(-(w-2*k) * (w-2*k)/2/t) + f
    return(1/sqrt(2*pi*t*t*t) * (f + w * exp(-w*w/2/t)))
}

#
# Density at lower barrier
#
# t: time (vector)
# v: drift
# a: upper barrier
# w: relative position of X(0) = z, w = z/a
# eps: required precision
#
fs = function(t, v, a, w, eps=sqrt(.Machine$double.eps))
{
    g = 1/a/a * exp(-v*a*w - v*v*t/2)
    return(g * fsw(t/a/a, w, eps/g))
}

#
# Examples
#
# Fs(t=1:1000, v=-0.1, a=100, w=0.4) for the distribution within 1-1000 ms
# fs(t=1:1000, v=-0.1, a=100, w=0.4) for the density
# Fs(t=1:1000, v=-(-0.1), a=100, w=1-0.4) for the upper barrier
# sigma = 0.1
# Fs(t=1:1000, v=-0.1/sigma, a=100/sigma, w=0.6) same for non-unit SD
#
