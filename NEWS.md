# fddm v0.6-0 (Release date: 2022-X-X)


## Other Changes
* character variables in med_dec transformed to factors.

# fddm v0.5-1 (Release date: 2022-3-15)

## New Features
* Added new method to compute the PDF of the DDM

## Bug Fixes
* Fixed HTML deprecation in documentation

## Other Changes
* Changed the function parameters in `dfddm()` that determine how the method of calculation is chosen.




# fddm v0.4-1 (Release date: 2021-12-10)

## Bug Fixes
* Fixed 'pi' issue with `Rcpp`




# fddm v0.4-0 (Release date: 2021-7-27)

## New Features
* Added the function `pfddm()` to calculate the cumulative distribution function (CDF) of the DDM; it uses the same model parameters as the existing `dfddm()` function.

## Bug Fixes
* Fixed issue where the log optional parameter was ignored and the non-log result was produced.




# fddm v0.3-3 (Release date: 2021-4-4)

## Bug Fixes
* Fixed issue on Solaris involving floor function and the int data type.
* Revamped handling of model parameters whose values are out-of-bounds.




# fddm v0.3-1 (Release date: 2021-3-21)

## New Features
* Added new parameter `sigma` to denote the diffusion coefficient of the Wiener process that underlies the DDM.

## Bug Fixes
* Corrected the required number of terms in the small-time approximation methods.
* `k` now consistently refers to the number of _individual_ terms instead of sometimes the number of _pairs_ of terms.

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
