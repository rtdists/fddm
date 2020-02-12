
library("devtools")
load_all(recompile = TRUE)

document()


sqrt(.Machine$double.eps)

### preliminary testing
rt <- seq(1,3,by=0.5)
resp <- rep(0,3)
eps <- 0.000001
dfddm(rt, resp, 1, 0.4, 0, 0.5, 0, TRUE, "Foster", "2017", "small", eps)
resp
eps

library("microbenchmark")
summary(microbenchmark(no_log = dfddm(rt, resp, 1, 0.4, 0, 0.5, 0, FALSE, "Foster", "2017", "small", eps),
                       ya_log = dfddm(rt, resp, 1, 0.4, 0, 0.5, 0, TRUE,  "Foster", "2017", "small", eps),
                       times = 10000, unit = "us"))

usethis::use_package("RWiener", type = "Suggests")

### random stuff
usethis::use_build_ignore("development.R")
usethis::use_build_ignore("examples/")
usethis::use_build_ignore("docs/")

usethis::use_travis()
usethis::use_readme_rmd()
