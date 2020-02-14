
library("devtools")
load_all(recompile = TRUE)

document()


sqrt(.Machine$double.eps)

### preliminary testing
rt <- seq(1,3,by=0.5)
eps <- 0.000001
dfddm(rt, 0, 1, 0.4, 0, 0.5, 0, TRUE, "Foster", "2017", "small", eps)



### random stuff
usethis::use_package("RWiener", type = "Suggests")

usethis::use_build_ignore("development.R")
usethis::use_build_ignore("examples/")
usethis::use_build_ignore("docs/")

usethis::use_travis()
usethis::use_readme_rmd()
