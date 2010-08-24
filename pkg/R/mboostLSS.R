###
# Generic implementation for multiple component LSS-models

### (glm/gam/m/black)boostLSS functions

mboostLSS <- function(formula, data = list(), families = list(),
                      control = boost_control(), weights = NULL, ...){
    cl <- match.call()
    fit <- mboostLSS_fit(formula = formula, data = data, families = families,
                         control = control, weights = weights, ...,
                         fun = mboost)
    attr(fit, "call") <- cl
    return(fit)
}

glmboostLSS <- function(formula, data = list(), families = list(),
                        control = boost_control(), weights = NULL, ...){
    cl <- match.call()
    fit <- mboostLSS_fit(formula = formula, data = data, families = families,
                         control = control, weights = weights, ...,
                         fun = glmboost)
    attr(fit, "call") <- cl
    return(fit)
}

gamboostLSS <- function(formula, data = list(), families = list(),
                        control = boost_control(), weights = NULL, ...){
    cl <- match.call()
    fit <- mboostLSS_fit(formula = formula, data = data, families = families,
                         control = control, weights = weights, ...,
                         fun = gamboost)
    attr(fit, "call") <- cl
    return(fit)
}

blackboostLSS <- function(formula, data = list(), families = list(),
                          control = boost_control(), weights = NULL, ...){
    cl <- match.call()
    fit <- mboostLSS_fit(formula = formula, data = data, families = families,
                         control = control, weights = weights, ...,
                         fun = blackboost)
    attr(fit, "call") <- cl
    return(fit)
}

###
# Todo:
# allow to specify a list of formulas (in formula)
mboostLSS_fit <- function(formula, data = list(), families = list(),
                          control = boost_control(), weights = NULL,
                          fun = mboost, ...){

    cl <- match.call()

    if (length(families) == 0)
        stop(sQuote("families"), " not specified")

    if ("offset" %in% names(list(...)))
        stop("Do not use argument ", sQuote("offset"),
             ". Please specify offsets via families")
    ### Use mu in family to specify offset in mu-family etc.

    if (is.list(formula)){
        if (!all(names(formula) %in% names(families)) ||
            length(unique(names(formula))) != length(names(families)))
            stop(sQuote("formula"), " can be either a formula or a named list",
                 " of formulas with same names as ",  sQuote("families"), ".")
        response <- eval(as.expression(formula[[1]][[2]]), envir = data,
                         enclos = environment(formula[[1]]))
    } else {
        response <- eval(as.expression(formula[[2]]), envir = data,
                         enclos = environment(formula))
        tmp <- vector("list", length = length(families))
        names(tmp) <- names(families)
        for (i in 1:length(tmp))
            tmp[[i]] <- formula
        formula <- tmp
    }

    mstop <- control$mstop
    control$mstop <- 1

    if (is.list(mstop)) {
        if (!all(names(mstop) %in% names(families)) ||
            length(unique(names(mstop))) != length(names(families)))
            stop(sQuote("mstop"), " can be either a scalar or a named list",
                 " of mstop values with same names as ",  sQuote("families"), "in ",
                 sQuote("boost_control"))
        mstop <- mstop[names(families)] ## sort in order of families
        mstop <- unlist(mstop)
    } else {
        if(length(mstop) != 1)
            stop(sQuote("mstop"), " can be either a scalar or a named list",
                 " of mstop values with same names as ",  sQuote("families"), "in ",
                 sQuote("boost_control"))
        mstop <- rep(mstop, length(families))
        names(mstop) <- names(families)
    }

    nu <- control$nu

    if (is.list(nu)) {
        if (!all(names(nu) %in% names(families)) ||
            length(unique(names(nu))) != length(names(families)))
            stop(sQuote("nu"), " can be either a scalar or a named list",
                 " of nu values with same names as ",  sQuote("families"), "in ",
                 sQuote("boost_control"))
        nu <- nu[names(families)] ## sort in order of families
        nu <- unlist(nu)
    } else {
        if(length(nu) != 1)
            stop(sQuote("nu"), " can be either a scalar or a named list",
                 " of nu values with same names as ",  sQuote("families"), "in ",
                 sQuote("boost_control"))
        nu <- rep(nu, length(families))
        names(nu) <- names(families)
    }

    if (is.list(control$risk) || is.list(control$center) || is.list(control$trace))
        stop(sQuote("risk"),", ", sQuote("center"), " and ", sQuote("trace") ,
             " cannot be lists in ", sQuote("boost_control"))

    trace <- control$trace
    control$trace <- FALSE

    w <- weights
    if (is.null(weights)) weights <- rep.int(1, NROW(response))
    weights <- mboost:::rescale_weights(weights)

    fit <- vector("list", length = length(families))
    names(fit) <- names(families)

    mods <- 1:length(fit)

    offset <- vector("list", length = length(mods))
    names(offset) <- names(families)
    for (j in mods){
        offset[[j]] <- families[[j]]@offset(y = response, w = weights)
        for (k in mods){
            for (l in mods){
                if (!is.null(offset[[l]]))
                    assign(names(offset)[l], families[[l]]@response(offset[[l]]),
                           environment(families[[k]]@ngradient))
            }
        }
    }
    for (j in mods){
        ## update value of nuisance parameters in families
        for (k in mods[-j]){
            if (!is.null(fit[[k]]))
                assign(names(fit)[k], fitted(fit[[k]], type = "response"),
                       environment(families[[j]]@ngradient))
        }
        ## use appropriate nu for the model
        control$nu <- nu[[j]]
        ## <FIXME> Do we need to recompute ngradient?
        fit[[j]] <- do.call(fun, list(formula[[names(families)[[j]]]],
                                      data = data, family = families[[j]],
                                      control=control, weights = w,
                                      ...))
    }
    if (trace)
        do_trace(current = 1, mstart = 0,
                 mstop = max(mstop),
                 risk = fit[[length(fit)]]$risk())

    ### set up a function for iterating boosting steps
    iBoost <- function(niter) {
        start <- sapply(fit, mstop)
        mvals <- vector("list", length(niter))
        for (j in 1:length(niter)){
            mvals[[j]] <- rep(start[j] + niter[j], max(niter))
            if (niter[j] > 0)
                mvals[[j]][1:niter[j]] <- (start[j] + 1):(start[j] + niter[j])
        }

        for (i in 1:max(niter)){
            for (j in mods){
                ## update value of nuisance parameters
                for (k in mods[-j])
                    assign(names(fit)[k], fitted(fit[[k]], type = "response"),
                           environment(get("ngradient", environment(fit[[j]]$subset))))
                ## update value of u, i.e. compute ngradient with new nuisance parameters
                evalq(u <- ngradient(y, fit, weights), environment(fit[[j]]$subset))
                ## update j-th component to "m-th" boosting step
                fit[[j]][mvals[[j]][i]]
            }
            if (trace){
                ## which is the current risk? rev() needed to get the last
                ## list element with maximum length
                whichRisk <- names(which.max(rev(lapply(lapply(fit, function(x) x$risk()), length))))
                do_trace(current = max(sapply(mvals, function(x) x[i])),
                         mstart = ifelse(firstRun, 0, max(start)),
                         mstop = ifelse(firstRun, max(niter) + 1, max(niter)),
                         risk = fit[[whichRisk]]$risk())
            }
        }
        return(TRUE)
    }

    #if (all(mstop == 1)){
    #    class(fit) <- c(paste(cl$fun, "LSS", sep=""), "mboostLSS")
    #    return(fit)
    #}

    if (any(mstop > 1)){
        ## actually go for initial mstop iterations!
        firstRun <- TRUE
        tmp <- iBoost(mstop - 1)
    }

    firstRun <- FALSE
    class(fit) <- c(paste(cl$fun, "LSS", sep=""), "mboostLSS")

    ### update to a new number of boosting iterations mstop
    ### i <= mstop means less iterations than current
    ### i >  mstop needs additional computations
    ### updates take place in THIS ENVIRONMENT,
    ### some models are CHANGED!
    attr(fit, "subset") <- function(i) {
        if (length(i) == 1)
            i <- rep(i, length(fit))

        msf <- mstop(fit)
        niter <- i - msf
        minStart <- min(msf[niter != 0], i[niter != 0])

        ## check if minStart bigger than mstop of parameters that are not touched
        #if (length(msf[niter == 0]) > 0 && minStart < min(msf[niter == 0]))
           #minStart <- min(msf[niter == 0])

        ## reduce models first (when necessary)
        if (any(msf > minStart)){
            cf <- class(fit)
            class(fit) <- "list" ## needed to use [] operator for lists

            #cat("processed parameters: ", paste(names(fit[msf > minStart]),
            #                                    collapse = ", "), "\n")

            lapply(fit[msf > minStart],
                   function(obj) obj$subset(minStart))

            ## remove additional boosting iterations from environments
            lapply(fit[msf > minStart], function(obj){
                   evalq({xselect <- xselect[1:mstop];
                          mrisk <- mrisk[1:mstop];
                          ens <- ens[1:mstop];
                          nuisance <- nuisance[1:mstop]},
                         environment(obj$subset))
               })

            class(fit) <- cf

            cat("Model first reduced to mstop = ", minStart, ".\n",
                "Now continue ...\n", sep ="")
        }

        ## now increase models (when necessary)
        if (any(i > minStart)){
            ## set negative values to 0
            ## (only applicable if some parameters do not need to be touched
            inc <- ifelse(i - minStart > 0, i - minStart, 0)
            tmp <- iBoost(inc)
        }

        mstop <<- i
    }
    return(fit)
}
