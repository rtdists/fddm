
data("med_dec", package = "fddm")
str(med_dec)

## number of participants per expertise condition:
aggregate(id ~ group, med_dec, function(x) length(unique(x)))

## number of trials per participant
aggregate(rt ~ group + id, med_dec, length)
