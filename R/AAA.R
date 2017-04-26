.onLoad <- function(libname, pkgname) {
    ### stabilize negative gradient by using a multiplicative factor based on
    ### the variance of the negative gradient
    options(gamboostLSS_stab_ngrad = FALSE)
}

# get rid of NOTEs in R CMD check for "undefined global functions or variables"
globalVariables(c("ngradient", "y", "combined_risk"))
