###
# Generic implementation for multiple component LSS-models

### (glm/gam/m/black)boostLSS functions

mboostLSS <- function(...)
    mboostLSS_fit(..., fun = mboost)

glmboostLSS <- function(...)
    mboostLSS_fit(..., fun = glmboost)

gamboostLSS <- function(...)
    mboostLSS_fit(..., fun = gamboost)

blackboostLSS <- function(...)
    mboostLSS_fit(..., fun = blackboost)

###
# Todo:
# allow to specify a list of formulas (in formula)
mboostLSS_fit <- function(formula, data, families = list(), ctrl = boost_control(),
                          weights = NULL, fun = mboost, ...){

    cl <- match.call()

    if ("offset" %in% names(list(...)))
        stop("Do not use argument ", sQuote("offset"),
        ". Please specify offsets via families")
    ### Use mu in family to specify offset in mu-family etc.

    mstop <- ctrl$mstop
    mstart <- 0
    ctrl$mstop <- 1

    trace <- ctrl$trace
    ctrl$trace <- FALSE

    w <<- eval(weights)

    fit <- vector("list", length = length(families))
    names(fit) <- names(families)

    mods <- 1:length(fit)
    for (j in mods){
        ## update value of nuisance parameters in families
        for (k in mods[-j]){
            if (!is.null(fit[[k]]))
                assign(names(fit)[k], fitted(fit[[k]], type = "response"),
                       environment(families[[j]]@ngradient))
        }
        ## <FIXME> Do we need to recompute ngradient?

        fit[[j]] <- fun(formula, data = data, family = families[[j]],
                        control=ctrl, weights = w, ...)
        fit2 <- fun(formula, data = data,
                    family = StudentTSigma(off = list(mu = fitted(fit[[1]]), df = 1)),
                    control=ctrl, weights = w, ...)
    }
    if (trace)
        mboost:::do_trace(1, mstop = mstart, risk = fit[[length(fit)]]$risk(),
                          width = mstop)
    if (mstop == 1){
        class(fit) <- c(paste(cl$fun, "LSS", sep=""), "mboostLSS")
        return(fit)
    }

    ## Loop
    for(m in (mstart + 2):(mstart + mstop)){
        for (j in mods){
            ## update value of nuisance parameters
            for (k in mods[-j])
                assign(names(fit)[k], fitted(fit[[k]], type = "response"),
                       environment(get("ngradient", environment(fit[[j]]$subset))))
            ## update value of u, i.e. compute ngradient with new nuisance parameters
            evalq(u <- ngradient(y, fit, weights), environment(fit[[j]]$subset))
            ## update j-th component to m boosting steps
            fit[[j]][m]
        }
        if (trace)
            mboost:::do_trace(m, mstop = mstart, risk = fit[[length(fit)]]$risk(),
                              width = mstop)
    }
    class(fit) <- c(paste(cl$fun, "LSS", sep=""), "mboostLSS")
    return(fit)
}
