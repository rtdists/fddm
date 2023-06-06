# prepare data
data(med_dec, package = "fddm")
med_dec <- med_dec[which(med_dec[["rt"]] >= 0), ] ## only use valid RTs
## select data from one participant
p1 <- med_dec[med_dec[["id"]] == 2 & med_dec[["group"]] == "experienced", ]
head(p1)


##----------------------------------------------
##        Easiest: Fitting using emmeans       -
##----------------------------------------------

## Because we use an ANOVA approach, we set orthogonal sum-to-zero contrasts
op <- options(contrasts = c('contr.sum', 'contr.poly'))

fit0 <- ddm(rt + response ~ classification*difficulty, data = p1)
summary(fit0)

## for more tests, we use emmeans:
if (requireNamespace("emmeans")) {
# for ANOVA table:
emmeans::joint_tests(fit0)
# get conditional main effects of difficulty (for each level of classification):
emmeans::joint_tests(fit0, by = "classification")

# get mean drift rates per condition:
em1 <- emmeans::emmeans(fit0, "difficulty", by = "classification")
em1
# compare mean drift rates per condition
pairs(em1)
update(pairs(em1), by = NULL, adjust = "holm")
}

options(op) # reset contrasts


##----------------------------------------------------------------
##              Fitting with custom parametrisation              -
##----------------------------------------------------------------

## one drift rate per classification by difficulty design cell
fit1 <- ddm(rt + response ~ 0 + classification:difficulty, data = p1)
summary(fit1)

## set default contrasts (just in case contrasts have been changed)
op <- options(contrasts = c('contr.treatment', 'contr.poly'))
## one drift rate "intercept" per classification condition (blast vs. non-blast)
## corresponding to first level of difficulty factor (easy)
## plus one further coefficient per classification condition corresponding to
## difference from "intercept" (hard - easy)
fit1b <- ddm(rt + response ~ 0 + classification + classification:difficulty,
             data = p1, args_optim = list(control = list(iter.max = 1000,
                                                         eval.max = 1000)))
summary(fit1b)
options(op) # reset contrasts

## set orthogonal sum-to-zero contrasts
op <- options(contrasts = c('contr.sum', 'contr.poly'))
## one drift rate "intercept" per classification condition (blast vs. non-blast)
## corresponding to mean drift rate for the classification condition
## plus one further coefficient per classification condition corresponding to
## difference from "intercept" (hard/easy - mean drift rate)
fit1c <- ddm(rt + response ~ 0 + classification + classification:difficulty,
             data = p1)
summary(fit1c)
options(op) ## reset contrasts

## all variants produce same fit, only meaning of parameters differs
logLik(fit1)
logLik(fit1b)
logLik(fit1c)
logLik(fit0) ## also model above

## all models estimate same drift rates, but in different parametrisation:
coef(fit1)  ## drift rates per design cell
## same drift rates based on fit1b:
c(coef(fit1b)[1:2],
  coef(fit1b)[1] + coef(fit1b)[3], coef(fit1b)[2] + coef(fit1b)[4])
## same drift rates based on fit1c:
c(coef(fit1c)[1] + coef(fit1c)[3], coef(fit1c)[2] + coef(fit1c)[4],
  coef(fit1c)[1] - coef(fit1c)[3], coef(fit1c)[2] - coef(fit1c)[4])

# we can estimate a model that freely estimates response bias
# (instead of fixing it at 0.5)
fit2 <- ddm(rt + response ~ 0 + classification:difficulty, bias = ~1, data = p1)
summary(fit2)
## Note: estimating bias only makes sense in situations like here where the
## response boundaries do not correspond to correct/incorrect but to the
## actual responses participants gave (here: blast vs. non-blast classification)

## now let's perform a likelihood ratio test to check if estimating response
## bias freely leads to a significant increase in model fit?
if (requireNamespace("lmtest")) { ## requires package lmtest
  lmtest::lrtest(fit1, fit2)
  ## does not look like it (p = 0.1691)
}


# we can also make a DDM parameter, such as boundary, depend on a numeric
# variable, such as block number
fit3 <- ddm(rt + response ~ 0 + classification:difficulty,
            boundary = ~ block, data = p1)
summary(fit3)

## does making boundary depend on block leads to a significant increase in model
## fit?
if (requireNamespace("lmtest")) { ## requires package lmtest
  lmtest::lrtest(fit1, fit3)
  ## does not look like it (p = 0.198)
}


##----------------------------------------------------------------
##              Fitting with optimization arguments              -
##----------------------------------------------------------------
## example of how to use your own initial values and bounds for the optimization
## of the coefficients (determined by the model matrices/formulas)
options(op) # reset contrasts

# this uses the default generated initial values and bounds
fitex0 <- ddm(rt + response ~ 0 + classification:difficulty, data = p1)

# we can see the number of coefficients (and thus the number of initial values
# and bounds) required if we wish to use our own
fitex0$optim_info$args_optim$init
fitex0$optim_info$args_optim$lo_bds
fitex0$optim_info$args_optim$up_bds
# -1.9270279  2.5408864 -1.9270279  0.3181987  1.7907438  0.1264035
# -Inf        -Inf      -Inf        -Inf       0          0
# Inf         Inf       Inf         Inf        Inf        0.461
# the first four coefficients are for the drift rate (given by the formula we
# provided), and the last two are for the boundary separation and non-decision
# time, respectively (note that the default is to estimate these two model
# parameters with a single coefficient).

# to use our own initial values, we can include them in the `args_optim` list,
# but note that they must be within the generated bounds
fitex1 <- ddm(rt + response ~ 0 + classification:difficulty, data = p1,
              args_optim = list(init = c(-1.5, 1, -1, 1, 1.5, 0.3)))

# to use our own bounds, we again can include them in the `args_optim` list,
# but note that they must contain the generated initial values
fitex2 <- ddm(rt + response ~ 0 + classification:difficulty, data = p1,
              args_optim = list(lo_bds = c(-20, -20, -20, -20, 0, 0),
                                up_bds = c(20, 20, 20, 20, 10, 0.5)))

# to use both our own initial values and bounds, we include them in the
# `args_optim` list, and the bounds must contain the initial values
fitex3 <- ddm(rt + response ~ 0 + classification:difficulty, data = p1,
              args_optim = list(init = c(-1.5, 1, -1, 1, 1.5, 0.3),
                                lo_bds = c(-20, -20, -20, -20, 0, 0),
                                up_bds = c(20, 20, 20, 20, 10, 0.5)))

# the only other option in the `args_optim` list is the control parameters
# directly used by the optimization function (e.g., the maximum number of
# iterations); here we'll set the maximum number of iterations and function
# evaluations
fitex4 <- ddm(rt + response ~ 0 + classification:difficulty, data = p1,
              args_optim = list(control = list(iter.max = 1000,
                                               eval.max = 1000)))
