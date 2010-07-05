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

    if (is.list(formula)){
        if (!all(names(formula) %in% names(families)) ||
            length(unique(names(formula))) != length(names(families)))
            stop(sQuote("formula"), " can be either a formula or a named list",
                 " of formulas with same names as ",  sQuote("families"), ".")
    } else {
        tmp <- vector("list", length = length(families))
        names(tmp) <- names(families)
        for (i in 1:length(tmp))
            tmp[[i]] <- formula
        formula <- tmp
    }


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
        fit[[j]] <- fun(formula[[names(families)[[j]]]], data = data, family = families[[j]],
                        control=ctrl, weights = w, ...)
    }
    if (trace)
        mboost:::do_trace(1, mstop = mstart, risk = fit[[length(fit)]]$risk(),
                          width = mstop)

    ### set up a function for iterating boosting steps
    iBoost <- function(niter) {
        start <- ifelse(mstart == 0, mstart + 2, mstart + 1)
        for(m in start:(mstart + niter)){
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
                                  width = niter)
        }
        mstart <<- mstart + niter
        return(TRUE)
    }

    if (mstop == 1){
        class(fit) <- c(paste(cl$fun, "LSS", sep=""), "mboostLSS")
        return(fit)
    }

    ### actually go for initial mstop iterations!
    tmp <- iBoost(mstop)

    class(fit) <- c(paste(cl$fun, "LSS", sep=""), "mboostLSS")

    ### update to a new number of boosting iterations mstop
    ### i <= mstop means less iterations than current
    ### i >  mstop needs additional computations
    ### updates take place in THIS ENVIRONMENT,
    ### some models are CHANGED!
    attr(fit, "subset") <- function(i) {
        if (i <= mstart){
            #stop("Use ", sQuote("model[i]"), " to decrease mstop")
            lapply(fit, function(a) a$subset(i))
            mstop <- i
        } else {
            tmp <- iBoost(i - mstart)
            mstop <- i
        }
    }
    return(fit)
}
