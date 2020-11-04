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
