# prepare data
data(med_dec, package = "fddm")
med_dec <- med_dec[which(med_dec[["rt"]] >= 0), ]
p1 <- med_dec[med_dec[["id"]] == 2 &
                      med_dec[["group"]] == "experienced", ]
p1[["response"]] <- as.factor(p1[["response"]])

# example fit with gradient information
ddm(rt + response ~ 0 + classification:difficulty, data = p1)

# example fit without gradient without gradient information
ddm(rt + response ~ 0 + classification:difficulty, data = p1,
    use_gradient = FALSE)
