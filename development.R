########## Do the following for every package update!!!
devtools::document()
devtools::build_vignettes()
devtools::build_manual()
# knit README.Rmd
devtools::build()
devtools::install()
# run R CMD check
devtools::check()


### TODO ###
--determine sl_thresh for second partials
--tests for hessian
update vignettes
--update other documentation
update example in README to include ddm() function
--reorder arguments pdf and derivs, v before a
--v1.0
dfddm docs link to fitting function
--add hessian slot (with vcov)
documentation for ddm function with list of slots (and describe a little bit)
maybe write "how to get started beginner's guide to fitting the DDM"
lingering issues:
  - standard errors are all the same (model matrix input is 1's and -1's so the main diagonal of the Hessian is always the same)

###
library("devtools")
devtools::load_all(recompile = TRUE)
devtools::load_all()
library("fddm")
devtools::document()

###################
# because Solaris and rtdists don't mix well
rhub::validate_email("singmann@gmail.com")
rhub::check_on_solaris(env_vars = c(`_R_CHECK_FORCE_SUGGESTS_` = "false"))
rhub::check_for_cran(path = "../fddm_0.4-1.tar.gz",email = "singmann@gmail.com")

sqrt(.Machine$double.eps)


### random stuff
use_testthat()
devtools::test()
devtools::check(vignettes = FALSE)
devtools::clean_dll()

usethis::use_package("ggnewscale", type = "Suggests")
usethis::use_rcpp_eigen()

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
