
<!-- README.md is generated from README.Rmd. Please edit that file -->

# fddm

<!-- badges: start -->

[![R build
status](https://github.com/rtdists/fddm/workflows/R-CMD-check/badge.svg)](https://github.com/rtdists/fddm/actions)
[![R-CMD-check](https://github.com/rtdists/fddm/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/rtdists/fddm/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

`fddm` provides the function `dfddm()`, which evaluates the density
function (or probability density function, PDF) for the Ratcliff
diffusion decision model (DDM) using different methods for approximating
the full PDF, which contains an infinite sum. `fddm` also provides the
family of functions `d*_dfddm()`, which evaluate the first-order partial
derivatives of the DDM density function with respect to the parameter
indicated by the `*` in the function name; the available parameters are
listed below. Similarly, `fddm` provides the family of functions
`d*2_dfddm()`, which evaluate the second-order partial derivatives of
the DDM density function for the same parameters. Based on the density
function and its partial derivatives, `fddm` provides the function
`ddm()`, which fits the DDM to provided data. `fddm` also provides the
function `pfddm()`, which evaluates the distribution function (or
cumulative distribution function, CDF) for the DDM using two different
methods for approximating the CDF.

Our implementation of the DDM has the following parameters: *a ϵ (0,
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
#>  $ id            : Factor w/ 37 levels "1","2","3","4",..: 2 2 2 2 2 2 2 2 2 2 ...
#>  $ group         : Factor w/ 3 levels "experienced",..: 1 1 1 1 1 1 1 1 1 1 ...
#>  $ block         : int  3 3 3 3 3 3 3 3 3 3 ...
#>  $ trial         : int  1 2 3 4 5 6 7 8 9 10 ...
#>  $ classification: Factor w/ 2 levels "blast","non-blast": 1 2 2 2 1 1 1 1 2 1 ...
#>  $ difficulty    : Factor w/ 2 levels "easy","hard": 1 1 2 2 1 1 2 2 1 2 ...
#>  $ response      : Factor w/ 2 levels "blast","non-blast": 1 2 1 2 1 1 1 1 2 1 ...
#>  $ rt            : num  0.853 0.575 1.136 0.875 0.748 ...
#>  $ stimulus      : Factor w/ 312 levels "blastEasy/AuerRod.jpg",..: 7 167 246 273 46 31 132 98 217 85 ...
```

### Easy Fitting with Built-in `ddm()`

The `ddm()` function fits the 5-parameter DDM to the user-supplied data
via maximum likelihood estimation. Each DDM parameter can be modeled
using R’s formula interface; the model parameters can either be fixed or
estimated, except for the drift rate which is always estimated.

We will demonstrate a simple example of how to fit the DDM to the `onep`
dataset from the above code chunks.

Because we use an ANOVA approach, we set orthogonal sum-to-zero
contrasts.

``` r
op <- options(contrasts = c('contr.sum', 'contr.poly'))
```

Now we can use the `ddm()` function to fit the DDM to the data. Note
that we are using formula notation to indicate the interaction between
variables for the drift rate. The first argument of the `ddm()` function
is the formula indicating how the drift rate should be modeled. By
default, the boundary separation and non-decision time are estimated,
and the initial bias and inter-trial variability in the drift rate are
held fixed.

``` r
fit0 <- ddm(rt + response ~ classification*difficulty, data = onep)
summary(fit0)
#> 
#> Call:
#> ddm(drift = rt + response ~ classification * difficulty, data = onep)
#> 
#> DDM fit with 3 estimated and 2 fixed distributional parameters.
#> Fixed: bias = 0.5, sv = 0 
#> 
#> drift coefficients (identity link):
#>                             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)                  -0.5924     0.1168  -5.073 3.91e-07 ***
#> classification1              -2.6447     0.1168 -22.647  < 2e-16 ***
#> difficulty1                   0.2890     0.1168   2.475   0.0133 *  
#> classification1:difficulty1  -1.4987     0.1168 -12.834  < 2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> boundary coefficients (identity link):
#>             Estimate Std. Error
#> (Intercept)    2.064      0.058
#> 
#> ndt coefficients (identity link):
#>             Estimate Std. Error
#> (Intercept)   0.3938      0.007
```

This output first shows the input to the function call and which
parameters are held fixed; this information is useful to verify that the
formula inputs to the `ddm()` function call were correct. For the model
of the drift rate that we input, we can see the estimates and summary
statistics for each coefficient. Below this, we can see the simple
estimates for the boundary separation and non-decision time (default
behavior).

We can reset the contrasts after fitting.

``` r
options(op) # reset contrasts
```

### Alternative Fitting Method with `nlminb()`

Although we strongly recommend using the `ddm()` function for fitting
the DDM to data because it is faster and more convenient, we will also
show how to use the probability density function in a manual
optimization setup. We further prepare the data by defining upper and
lower responses and the correct response bounds.

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
              control = list(iter.max = 300, eval.max = 300),
              rt = onep[["rt"]], resp = onep[["resp"]], truth = onep[["truth"]],
              # limits:   vu,   vl,   a,                t0, w,  sv
              lower = c(-Inf, -Inf, .01,                 0, 0,   0),
              upper = c( Inf,  Inf, Inf, min(onep[["rt"]]), 1, Inf))
fit
#> $par
#> [1]  5.6813148 -2.1886625  2.7909164  0.3764463  0.4010114  2.2812999
#> 
#> $objective
#> [1] 42.47181
#> 
#> $convergence
#> [1] 0
#> 
#> $iterations
#> [1] 242
#> 
#> $evaluations
#> function gradient 
#>      266     1723 
#> 
#> $message
#> [1] "relative convergence (4)"
```
