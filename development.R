
library("devtools")
devtools::load_all(recompile = TRUE)
library("fddm")
devtools::document()

###################


sqrt(.Machine$double.eps)

### preliminary testing
rt <- seq(1,3,by=0.5)
eps <- 1e-6
dfddm(rt, 1, 1, 0.4, 0, 0.5, 0, T, "Gondan", "2017", "small", eps)
dfddm(rt = numeric(), response = "lower", a = 1, v = -2, t0 = 0.3)



### random stuff
use_testthat()
devtools::test()

usethis::use_package("ggplot2", type = "Suggests")

usethis::use_build_ignore("development.R")
usethis::use_build_ignore("examples/")
usethis::use_build_ignore("docs/")

usethis::use_travis()
usethis::use_readme_rmd()
usethis::use_github_action_check_standard()

use_vignette("Validity")
devtools::build_vignettes()
devtools::install(build_vignettes = TRUE)
# html_document aesthetic options that I can't use........
output:
  html_document:
    theme: paper
    highlight: default
    keep_md: true
    toc: false
    code_folding: hide
    fig_width: 16
    fig_height: 9
    fig_caption: true

med_dec <- read.csv("inst/extdata/medical_dm.csv", stringsAsFactors = FALSE)
save(med_dec, file = "data/med_dec.rda", compress = "bzip2", compression_level = 9)
tools::checkRdaFiles("data/med_dec.rda")
usethis::use_data(med_dec, overwrite = TRUE, version = 3, name = "med_dec",
                  compress = "bzip2", compression_level = 9)
