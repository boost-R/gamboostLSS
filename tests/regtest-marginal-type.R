require("gamboostLSS")

### check detection of get_marginal_type for continuous
stopifnot(gamboostLSS:::get_marginal_type(GaussianLSS()) == 
            "continuous")
stopifnot(gamboostLSS:::get_marginal_type(GammaLSS()) == 
            "continuous")
stopifnot(gamboostLSS:::get_marginal_type(LogNormalLSS()) == 
            "continuous")
stopifnot(gamboostLSS:::get_marginal_type(BetaLSS()) == 
            "continuous")

### check detection of get_marginal_type for discrete
stopifnot(gamboostLSS:::get_marginal_type(NBinomialLSS()) == "discrete")
stopifnot(gamboostLSS:::get_marginal_type(ZIPoLSS()) == "discrete")
stopifnot(gamboostLSS:::get_marginal_type(ZINBLSS()) == "discrete")

### get_marginal_type throws errors correctly
tryCatch(gamboostLSS:::get_marginal_type(list()),
         error = function(e) cat("Unknown family throws error: OK\n"))

### check detection of get_marginal_case
stopifnot(gamboostLSS:::get_marginal_case(GaussianLSS(), GaussianLSS()) == 
            "continuous_continuous")
stopifnot(gamboostLSS:::get_marginal_case(GaussianLSS(), ZIPoLSS()) == 
            "continuous_discrete")
stopifnot(gamboostLSS:::get_marginal_case(ZIPoLSS(), GaussianLSS()) == 
            "discrete_continuous")