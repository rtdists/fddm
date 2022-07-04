# prepare data
data(med_dec, package = "fddm")
med_dec <- med_dec[which(med_dec[["rt"]] >= 0), ] ## only use valid RTs
## select data from one participant
p1 <- med_dec[med_dec[["id"]] == 2 & med_dec[["group"]] == "experienced", ]
head(p1)

## one drift rate per classification by difficulty design cell 
fit1 <- ddm(rt + response ~ 0 + classification:difficulty, data = p1)
summary(fit1)
fit1

## set default contrasts (just in case contrasts have been changed)
op <- options(contrasts=c('contr.treatment', 'contr.poly'))
## one drift rate "intercept" per classification condition (blast vs. non-blast) 
## corresponding to first level of difficulty factor (easy)
## plus one further coefficient per classification condition corresponding to 
## difference from "intercept" (hard - easy)
fit1b <- ddm(rt + response ~ 0 + classification + classification:difficulty, 
             data = p1)
fit1b
options(op) # reset contrasts

## set orthogonal sum-to-zero contrasts
op <- options(contrasts=c('contr.sum', 'contr.poly'))
## one drift rate "intercept" per classification condition (blast vs. non-blast) 
## corresponding to mean drift rate for the classification condition 
## plus one further coefficient per classification condition corresponding to 
## difference from "intercept" (hard/easy - mean drift rate)
fit1c <- ddm(rt + response ~ 0 + classification + classification:difficulty, 
             data = p1)
fit1c
options(op) ## reset contrasts

## all three variants produce same fit, only meaning of parameters differs 
logLik(fit1)
logLik(fit1b)
logLik(fit1c)

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
fit2
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
fit3

## does making boundary depend on block leads to a significant increase in model
## fit?
if (requireNamespace("lmtest")) { ## requires package lmtest
  lmtest::lrtest(fit1, fit3)
  ## does not look like it (p = 0.198)
}

