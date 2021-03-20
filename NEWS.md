# fddm v0.3-1 (Release date: 2021-3-21)

## New Features
* Added new parameter `sigma` to denote the diffusion coefficient of the Wiener process that underlies the DDM.

## Bug Fixes
* Amended the required number of terms and summation functions.

## Other Changes
* Changed default value of `max_terms_large` to 1. This affects the default method (`n_terms_small = "SWSE", scale = "both"`).
* Revamped the error handling to make it less intrusive.
* Removed function call in `R`; now directly calls `C++`.
* Updated vignettes.




# fddm v0.2-2 (Release date: 2020-11-03)

## Other Changes
* Added warning about old 32-bit compilers exhibiting abnormal behavior for the large-time approximation.




# fddm v0.2-1 (Release date: 2020-10-09)

## Bug Fixes
* The num_funcs.cpp functions are now prevented from possible overflow.
* Checks are now in place for empty strings as inputs to the parameters `n_terms_small`, `summation_small`, and `scale`.




# fddm v0.1-2 (Release date: 2020-10-07)

## New Features
* New default method, callable using the arguments `n_terms_small = "SWSE", scale = "both"`.

## Bug Fixes
* Density functions now return strictly non-negative values.

## Other Changes
* If input response times is an empty vector, then an empty vector is returned.
* Fine-tuned the summation functions used with the SWSE method
* Added URLs to DESCRIPTION




# fddm v0.1-1 (Release date: 2020-07-10)

Initial CRAN release.
