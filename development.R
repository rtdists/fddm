########## Do the following for every package update!!!
devtools::document()
devtools::build_vignettes()
devtools::build_manual()
# knit README.Rmd using code below:
devtools::build_readme()
devtools::build()
devtools::install()
# run R CMD check
devtools::check(cran = TRUE)


##### R CMD Check Warnings/Notes #####

### Warnings
#- ‘qpdf’ is needed for checks on size reduction of PDFs
#---> I looked this up, and we can safely ignore this

### Notes
#- Specified C++11: please drop specification unless essential
#---> it's essential for Rcpp to work properly

#- installed size is 15.8Mb
    # sub-directories of 1Mb or more:
    #   doc    1.4Mb
    #   libs  13.5Mb
#---> dunno what we can do about that

#- Namespace in Imports field not imported from: ‘RcppEigen’
    # All declared Imports should be used.
#---> it's ok; we use RcppEigen to do matrix multiplication, but it's an
      # operator (*) instead of a function() so it looks like we don't use it
#-------------------------------------

### TODO ###

### Maybe Later ###
# rerun benchmark vignette and verify the comments on the plots
# example vignette- add fitting entire dataset with ddm() and some simple analysis
# write vignette "how to get started beginner's guide to fitting the DDM"
# there is also the possibility of adding more link functions (instead of just the identity)


###
library("devtools")
devtools::load_all(recompile = TRUE)
devtools::load_all()
library("fddm")
devtools::document()

###################
# because Solaris and rtdists don't mix well
rhub::rhub_doctor()
rhub::rhub_check() 

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
usethis::use_github_action("check-standard")

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
