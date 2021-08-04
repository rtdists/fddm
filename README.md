
<!-- README.md is generated from README.Rmd. Please edit that file -->

# fddm

<!-- badges: start -->

[![Travis build
status](https://travis-ci.org/rtdists/fddm.svg?branch=master)](https://travis-ci.org/rtdists/fddm)
[![R build
status](https://github.com/rtdists/fddm/workflows/R-CMD-check/badge.svg)](https://github.com/rtdists/fddm/actions)
<!-- badges: end -->

`fddm` provides function `dfddm()` which evaluates the density function
(or probability density function, PDF) for the Ratcliff diffusion
decision model (DDM) using different methods for approximating the full
PDF, which contains an infinite sum. Our implementation of the DDM has
the following parameters: *a ϵ (0,
<font style="vertical-align: middle;" size="5em">∞</font>)* (threshold
separation), *v ϵ
(-<font style="vertical-align: middle;" size="5em">∞</font>,
<font style="vertical-align: middle;" size="5em">∞</font>)* (drift
rate), *t<sub>0</sub> ϵ \[0,
<font style="vertical-align: middle;" size="5em">∞</font>)*
(non-decision time/response time constant), *w ϵ (0, 1)* (relative
starting point), *sv ϵ (0,
<font style="vertical-align: middle;" size="5em">∞</font>)*
(inter-trial-variability of drift), and *sigma ϵ (0,
<font style="vertical-align: middle;" size="5em">∞</font>)* (diffusion
coefficient of the underlying Wiener Process).

## Installation

You can install the released version of fddm from
[CRAN](https://CRAN.R-project.org) with:

    install.packages("fddm")

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("rtdists/fddm")
```

## Example

As a preliminary example, we will fit the DDM to the data from one
participant in the `med_dec` data that comes with `fddm`. This dataset
contains the accuracy condition reported in Trueblood et al. (2018),
which investigates medical decision making among medical professionals
(pathologists) and novices (i.e., undergraduate students). The task of
participants was to judge whether pictures of blood cells show cancerous
cells (i.e., blast cells) or non-cancerous cells (i.e., non-blast
cells). The dataset contains 200 decisions per participant, based on
pictures of 100 true cancerous cells and pictures of 100 true
non-cancerous cells. Here we use the data collected from the trials of
one experienced medical professional (pathologist). First, we load the
`fddm` package, remove any invalid responses from the data, and select
the individual whose data we will use for fitting.

``` r
library("fddm")
data(med_dec, package = "fddm")
med_dec <- med_dec[which(med_dec[["rt"]] >= 0), ]
onep <- med_dec[ med_dec[["id"]] == "2" & med_dec[["group"]] == "experienced", ]
str(onep)
#> 'data.frame':    200 obs. of  9 variables:
#>  $ id            : int  2 2 2 2 2 2 2 2 2 2 ...
#>  $ group         : chr  "experienced" "experienced" "experienced" "experienced" ...
#>  $ block         : int  3 3 3 3 3 3 3 3 3 3 ...
#>  $ trial         : int  1 2 3 4 5 6 7 8 9 10 ...
#>  $ classification: chr  "blast" "non-blast" "non-blast" "non-blast" ...
#>  $ difficulty    : chr  "easy" "easy" "hard" "hard" ...
#>  $ response      : chr  "blast" "non-blast" "blast" "non-blast" ...
#>  $ rt            : num  0.853 0.575 1.136 0.875 0.748 ...
#>  $ stimulus      : chr  "blastEasy/BL_10166384.jpg" "nonBlastEasy/16258001115A_069.jpg" "nonBlastHard/BL_11504083.jpg" "nonBlastHard/MY_9455143.jpg" ...
```

We further prepare the data by defining upper and lower responses and
the correct response bounds.

``` r
onep[["resp"]] <- ifelse(onep[["response"]] == "blast", "upper", "lower")
onep[["truth"]] <- ifelse(onep[["classification"]] == "blast", "upper", "lower")
```

For fitting, we need a simple likelihood function; here we will use a
straightforward log of sum of densities of the study responses and
associated response times. This log-likelihood function will fit the
standard parameters in the DDM, but it will fit two versions of the
drift rate *v*: one for when the correct response is `"blast"`
(*v<sub>u</sub>*), and another for when the correct response is
`"non-blast"` (*v<sub>l</sub>*). A detailed explanation of the
log-likelihood function is provided in the Example Vignette
(`vignette("example", package = "fddm")`). Note that this likelihood
function returns the negative log-likelihood as we can simply minimize
this function to get the maximum likelihood estimate.

``` r
ll_fun <- function(pars, rt, resp, truth) {
  v <- numeric(length(rt))

  # the truth is "upper" so use vu
  v[truth == "upper"] <- pars[[1]]
  # the truth is "lower" so use vl
  v[truth == "lower"] <- pars[[2]]

  dens <- dfddm(rt = rt, response = resp, a = pars[[3]], v = v,
                t0 = pars[[4]], w = pars[[5]], sv = pars[[6]], log = TRUE)

  return( ifelse(any(!is.finite(dens)), 1e6, -sum(dens)) )
}
```

We then pass the data and log-likelihood function to an optimization
function with the necessary additional arguments. As we are using the
optimization function `nlminb` for this example, the first argument must
be the initial values of our DDM parameters that we want optimized.
These are input in the order: *v<sub>u</sub>*, *v<sub>l</sub>*, *a*,
*t<sub>0</sub>*, *w*, and *sv*; we also need to define upper and lower
bounds for each of the parameters. Fitting the DDM to this dataset is
basically instantaneous using this setup.

``` r
fit <- nlminb(c(0, 0, 1, 0, 0.5, 0), objective = ll_fun,
              rt = onep[["rt"]], resp = onep[["resp"]], truth = onep[["truth"]],
              # limits:   vu,   vl,   a,  t0, w,  sv
              lower = c(-Inf, -Inf, .01,   0, 0,   0),
              upper = c( Inf,  Inf, Inf, min(onep[["rt"]]), 1, Inf))
fit
#> $par
#> [1]  5.6813074 -2.1886617  2.7909132  0.3764465  0.4010115  2.2813001
#> 
#> $objective
#> [1] 42.47181
#> 
#> $convergence
#> [1] 0
#> 
#> $iterations
#> [1] 55
#> 
#> $evaluations
#> function gradient 
#>       82      389 
#> 
#> $message
#> [1] "relative convergence (4)"
```
