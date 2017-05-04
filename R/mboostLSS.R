###
# Generic implementation for multiple component LSS-models

### (glm/gam/m/black)boostLSS functions

mboostLSS <- function(formula, data = list(), families = GaussianLSS(),
                      control = boost_control(), weights = NULL,
                      method = c("cyclic", "noncyclic"), ...){

    cl <- match.call()
    if(is.null(cl$families))
        cl$families <- families
    method <- match.arg(method)

    fit <- mboostLSS_fit(formula = formula, data = data, families = families,
                         control = control, weights = weights, ...,
                         fun = mboost, funchar = "mboost", call = cl,
                         method = method)
    return(fit)
}

glmboostLSS <- function(formula, data = list(), families = GaussianLSS(),
                        control = boost_control(), weights = NULL,
                        method = c("cyclic", "noncyclic"), ...){

    cl <- match.call()
    if(is.null(cl$families))
        cl$families <- families
    method <- match.arg(method)

    fit <- mboostLSS_fit(formula = formula, data = data, families = families,
                         control = control, weights = weights, ...,
                         fun = glmboost, funchar = "glmboost", call = cl,
                         method = method)

    return(fit)
}

gamboostLSS <- function(formula, data = list(), families = GaussianLSS(),
                        control = boost_control(), weights = NULL,
                        method = c("cyclic", "noncyclic"), ...){

    cl <- match.call()
    if(is.null(cl$families))
        cl$families <- families
    method <- match.arg(method)

    fit <- mboostLSS_fit(formula = formula, data = data, families = families,
                         control = control, weights = weights, ...,
                         fun = gamboost, funchar = "gamboost", call = cl,
                         method = method)

    return(fit)
}

blackboostLSS <- function(formula, data = list(), families = GaussianLSS(),
                          control = boost_control(), weights = NULL,
                          method = c("cyclic", "noncyclic"), ...){
    cl <- match.call()
    if(is.null(cl$families))
        cl$families <- families
    method <- match.arg(method)

    fit <- mboostLSS_fit(formula = formula, data = data, families = families,
                         control = control, weights = weights, ...,
                         fun = blackboost, funchar = "blackboost", call = cl,
                         method = method)
    return(fit)
}


### work horse for fitting (glm/gam/m/black)boostLSS models
mboostLSS_fit <- function(formula, data = list(), families = GaussianLSS(),
                          control = boost_control(), weights = NULL,
                          fun = mboost, funchar = "mboost", call = NULL,
                          method, ...){

    if (length(families) == 0)
        stop(sQuote("families"), " not specified")

    if ("offset" %in% names(list(...)))
        stop("Do not use argument ", sQuote("offset"),
             ". Please specify offsets via families")
    ### Use mu in "families" to specify offset in mu-family etc.

    if (is.list(formula)){
        if (!all(names(formula) %in% names(families)) ||
            length(unique(names(formula))) != length(names(families)))
            stop(sQuote("formula"), " can be either a formula or a named list",
                 " of formulas with same names as ",  sQuote("families"), ".")
        ynames <- sapply(formula, function(fm) as.character(fm[[2]]))
        if (length(unique(ynames)) > 1)
            warning("responses differ for the components")
        response <- lapply(formula, function(fm)
            eval(as.expression(fm[[2]]), envir = data,
                 enclos = environment(fm)))
        #response <- eval(as.expression(formula[[1]][[2]]), envir = data,
        #                 enclos = environment(formula[[1]]))
    } else {
        response <- eval(as.expression(formula[[2]]), envir = data,
                         enclos = environment(formula))
        tmp <- vector("list", length = length(families))
        names(tmp) <- names(families)
        for (i in 1:length(tmp))
            tmp[[i]] <- formula
        formula <- tmp
    }

    mstop <- mstoparg <- control$mstop
    control$mstop <- 0
    if (method == "cyclic")
        mstop <- check(mstop, "mstop", names(families))


    nu <- control$nu
    nu <- check(nu, "nu", names(families))

    if (is.list(control$risk) || is.list(control$center) || is.list(control$trace))
        stop(sQuote("risk"),", ", sQuote("center"), " and ", sQuote("trace") ,
             " cannot be lists in ", sQuote("boost_control"))

    trace <- control$trace
    control$trace <- FALSE

    ## generate adequate model weights 
    w <- weights
    
    if (is.null(weights)){
      if (!is.list(response)) {
        weights <- rep.int(1, NROW(response))
        # expand weights if the response is a matrix (functional response)
        if(funchar == "FDboost" && !is.null(dim(response)) && !any(dim(response) == 1))
          weights <- rep.int(weights, ncol(response))
        
      } else {
        weights <- rep.int(1, NROW(response[[1]]))
        # expand weights if the response is a matrix (functional response)
        if(funchar == "FDboost" && !is.null(dim(response[[1]])) && !any(dim(response[[1]]) == 1))
          weights <- rep.int(weights, ncol(response[[1]]))
      }
    }
    
    weights <- rescale_weights(weights)
    
    ## set up timeformula for FDboost 
    if (funchar == "FDboost"){
      
      # get timeformula from dots 
      dots <- list(...)
      timeformula <- dots$timeformula
      dots$timeformula <- NULL
      
      # deal with argument timeformula in case of FDboost()
      # timeformula is named list with names according to 
      # distribution parameters of families
      timeformula <- check_timeformula(timeformula, families)
    }

    fit <- vector("list", length = length(families))
    names(fit) <- names(families)

    mods <- 1:length(fit)

    offset <- vector("list", length = length(mods))
    names(offset) <- names(families)
    for (j in mods){
        if (!is.list(response)) {
          response <- check_y_family(response, families[[j]], 
                                     allow_matrix = (funchar == "FDboost"))
          offset[[j]] <- families[[j]]@offset(y = if(funchar != "FDboost") response else c(response), 
                                              w = weights)
        } else {
          response[[j]] <- check_y_family(response[[j]], families[[j]], 
                                          allow_matrix = (funchar == "FDboost"))
          offset[[j]] <- families[[j]]@offset(y = if(funchar != "FDboost") response[[j]] else c(response[[j]]), 
                                              w = weights)
        }
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
            if (!is.null(fit[[k]])) ## use fitted.mboost() as fitted.FDboost() returns a matrix  
                assign(names(fit)[k], families[[k]]@response(fitted.mboost(fit[[k]])),
                       environment(families[[j]]@ngradient))
        }
        ## use appropriate nu for the model
        control$nu <- nu[[j]]
        ## <FIXME> Do we need to recompute ngradient?
        if(funchar != "FDboost"){
          fit[[j]] <- do.call(fun, list(formula[[names(families)[[j]]]],
                                        data = data, family = families[[j]],
                                        control = control, weights = w,
                                        ...))
        }else{
          fit[[j]] <- do.call(fun, c(list(formula[[names(families)[[j]]]],
                                        timeformula = timeformula[[names(families)[[j]]]],
                                        data = data, family = families[[j]],
                                        control = control, weights = w, 
                                        # always use scalar offset, as offsets are treated within the Family
                                        offset = "scalar"), 
                                     dots))
        }

    }

    iBoost <- function(niter, method) {

        start <- sapply(fit, mstop)
        # initialize combined_risk
        combined_risk <- NA

        if (method == "noncyclic") {
            ### noncyclical fitting ###
            # this is the case for boosting from the beginning
            if (is.null(attr(fit, "combined_risk")) | niter == 0) {
                combined_risk <- vapply(fit, risk, numeric(1))
            } else {
                combined_risk <- attr(fit, "combined_risk")()
            }

            best <- which(names(fit) == tail(names(combined_risk), 1))
        } else {
            ### cyclical fitting ###
            mvals <- vector("list", length(niter))
            for (j in 1:length(niter)) {
                mvals[[j]] <- rep(start[j] + niter[j], max(niter))
                if (niter[j] > 0)
                    mvals[[j]][1:niter[j]] <- (start[j] + 1):(start[j] + niter[j])
            }
        }

        ENV <- lapply(mods, function(j) environment(fit[[j]]$subset))
        ## main loop starts here ##
        for (i in 1:max(niter)){
            if (method == "noncyclic") {
                ### noncyclical fitting ###

                ## update value of nuisance parameters
                ## use response(fitted()) as this is much quicker than fitted(, type = response)
                for( k in mods[-best]) {
                  ## use fitted.mboost() as fitted.FDboost() returns a matrix
                    assign(names(fit)[best], families[[best]]@response(fitted.mboost(fit[[best]])),
                           environment(get("ngradient", environment(fit[[k]]$subset))))
                }

                risks <- numeric(length(fit))
                for(b  in 1:length(fit)){
                    st <- mstop(fit[[b]])
                    mstop(fit[[b]]) = st + 1
                    risks[b] <- tail(risk(fit[[b]]), 1)
                    #evalq({riskfct(y, fit, weights)}, envir = ENV[[b]])
                    #risks[b] <- ENV[[b]][["riskfct"]](ENV[[b]][["y"]], ENV[[b]][["fit"]], ENV[[b]][["weights"]])
                    fit[[b]][st]

                    ## fit[[b]][st] is not enough to reduce the model back to beginning, so
                    ## so all these values have to be reduced, so that they are calculated
                    ## correctly the next time
                    evalq({
                        xselect <- xselect[seq_len(mstop)];
                        mrisk <- mrisk[seq_len(mstop + 1)];
                        ens <- ens[seq_len(mstop)];
                        nuisance <- nuisance[seq_len(mstop)]
                    },
                    environment(fit[[b]]$subset))
                }

                best <- which.min(risks)

                ## update value of u, i.e. compute ngradient with new nuisance parameters
                evalq({u <- ngradient(y, fit, weights)}, ENV[[best]])
                #ENV[[best]][["u"]] <- ENV[[best]][["ngradient"]](ENV[[best]][["y"]], ENV[[best]][["fit"]], ENV[[best]][["weights"]])

                ## update selected component by 1
                fit[[best]][mstop(fit[[best]]) + 1]

                ## update risk list
                combined_risk[(length(combined_risk) + 1)] <- tail(risk(fit[[best]]), 1)
                names(combined_risk)[length(combined_risk)] <- names(fit)[best]
                combined_risk <<- combined_risk

            } else {
                ### cyclical fitting ###
                for (j in mods){
                    ## update value of nuisance parameters
                    ## use response(fitted()) as this is much quicker than fitted(, type = response)
                    for (k in mods[-j]) ## use fitted.mboost() as fitted.FDboost() returns a matrix
                        assign(names(fit)[k], families[[k]]@response(fitted.mboost(fit[[k]])),
                               environment(get("ngradient", environment(fit[[j]]$subset))))
                    ## update value of u, i.e. compute ngradient with new nuisance parameters

                    ENV[[j]][["u"]] <- ENV[[j]][["ngradient"]](ENV[[j]][["y"]], ENV[[j]][["fit"]],
                                                               ENV[[j]][["weights"]])
                    # same as:
                    # evalq(u <- ngradient(y, fit, weights), environment(fit[[j]]$subset))

                    ## update j-th component to "m-th" boosting step
                    fit[[j]][mvals[[j]][i]]
                }
            }

            if (trace){
                if (method == "noncyclic") {
                    do_trace(current = length(combined_risk) - length(fit), mstart = sum(start),
                             risk = combined_risk[-(1:length(fit))], mstop = niter)
                } else {
                    ## which is the current risk? rev() needed to get the last
                    ## list element with maximum length
                    whichRisk <- names(which.max(rev(lapply(fit, function(x) length(risk(x))))))
                    do_trace(current = max(sapply(mvals, function(x) x[i])),
                             mstart = max(start),
                             mstop = max(niter),
                             risk = risk(fit[[whichRisk]])[-1])
                }

            }
        }
        return(TRUE)
    }

    if (any(mstop > 0))
        tmp <- iBoost(mstop, method = method)

    class(fit) <- c(paste0(funchar, "LSS"), "mboostLSS")
    if(method != "cyclic"){
        class(fit) <- c("nc_mboostLSS", class(fit))
    }

    ### update to a new number of boosting iterations mstop
    ### i <= mstop means less iterations than current
    ### i >  mstop needs additional computations
    ### updates take place in THIS ENVIRONMENT,
    ### some models are CHANGED!
    if(method == "cyclic") {
        attr(fit, "subset") <- function(i) {

            i <- check(i, "mstop", names(families))

            msf <- mstop(fit)
            niter <- i - msf
            if (all(niter == 0)) {
                ## make nothing happen with the model
                minStart <- max(msf)
            } else {
                minStart <- min(msf[niter != 0], i[niter != 0])
            }

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
                    evalq({xselect <- xselect[seq_len(mstop)];
                    mrisk <- mrisk[seq_len(mstop + 1)];
                    ens <- ens[seq_len(mstop)];
                    nuisance <- nuisance[seq_len(mstop)]},
                    environment(obj$subset))
                })

                class(fit) <- cf
                if (trace)
                    cat("Model first reduced to mstop = ", minStart, ".\n",
                        "Now continue ...\n", sep ="")
            }

            ## now increase models (when necessary)
            if (any(i > minStart)){
                ## set negative values to 0
                ## (only applicable if some parameters do not need to be touched
                inc <- ifelse(i - minStart > 0, i - minStart, 0)
                tmp <- iBoost(inc, method = method)
            }

            mstop <<- i
        }
    }
    else {
        attr(fit, "subset") <- function(i) {
            msf <- sum(mstop(fit))
            niter <- i - msf

            ## check if minStart bigger than mstop of parameters that are not touched
            #if (length(msf[niter == 0]) > 0 && minStart < min(msf[niter == 0]))
            #minStart <- min(msf[niter == 0])

            ## reduce models first (when necessary)
            if (niter < 0 ){
                cf <- class(fit)
                d = length(fit)
                class(fit) <- "list" ## needed to use [] operator for lists
                #cat("processed parameters: ", paste(names(fit[msf > minStart]),
                #                                    collapse = ", "), "\n")
                #reduce the combined risk values
                combined_risk <- attr(fit, "combined_risk")()
                combined_risk <<- attr(fit, "combined_risk")()[seq_len(i + d)]
                new_stop_value <- table(names(attr(fit, "combined_risk")())) - 1


                for(o in names(new_stop_value)){
                    fit[[o]]$subset(new_stop_value[o])
                }

                ## remove additional boosting iterations from environments
                lapply(fit, function(obj){
                    evalq({xselect <- xselect[seq_len(mstop)];
                    mrisk <- mrisk[seq_len(mstop + 1)];
                    ens <- ens[seq_len(mstop)];
                    nuisance <- nuisance[seq_len(mstop)]},
                    environment(obj$subset))
                })

                #update ALL nuisance parameters to last update
                ENV <- lapply(mods, function(j) environment(fit[[j]]$subset))
                for(j in names(new_stop_value)){
                    for( k in setdiff(names(new_stop_value), j)){
                        assign(k, families[[k]]@response(fitted.mboost(fit[[k]])),
                               environment(get("ngradient", environment(fit[[j]]$subset))))
                    }
                }

                for(k in mods){
                    evalq(u <- get("ngradient")(get("y"), fit, weights), ENV[[k]])
                }


                class(fit) <- cf
            }

            ## now increase models (when necessary)
            else if (niter > 0){
                tmp <- iBoost(niter, method = method)
            }

            mstop <<- i
        }
    }

    ## make call in submodels nicer:
    cl <- call
    cl[[1]] <- as.name(gsub("LSS", "", cl[[1]]))
    names(cl)[names(cl) == "families"] <- "family"
    for (i in 1:length(fit)) {
        fit[[i]]$call <- cl
        ## <FIXME> This is not really nice
        fit[[i]]$call$family <- families[[i]]
    }

    attr(fit, "(weights)") <- weights  ## attach weights used for fitting

    ## update to new weights; just a fresh start
    attr(fit, "update") <- function(weights = NULL, oobweights = NULL,
                                    risk = NULL, trace = NULL, mstop = NULL) {
        if (is.null(mstop)) {
            control$mstop <- mstoparg
        } else {
            control$mstop <- mstop
        }
        if (!is.null(risk))
            control$risk <- risk
        if (!is.null(trace))
            control$trace <- trace
        ## re-use user specified offset only
        ## (since it depends on weights otherwise)
        ## this is achieved via a re-evaluation of the families argument
        
        if(funchar != "FDboost"){
          mboostLSS_fit(formula = formula, data = data,
                        families = eval(call[["families"]]), weights = weights,
                        control = control, fun = fun, funchar = funchar,
                        call = call, oobweights = oobweights,
                        method = method)
        }else{
          mboostLSS_fit(formula = formula, data = data,
                        families = eval(call[["families"]]), weights = weights,
                        control = control, fun = fun, funchar = funchar,
                        call = call, oobweights = oobweights,
                        method = method, 
                        timeformula = timeformula[[names(families)[[j]]]])
        }
        
    }
    attr(fit, "control") <- control
    attr(fit, "call") <- call
    attr(fit, "data") <- data
    attr(fit, "families") <- families
    if(method != "cyclic")
        attr(fit, "combined_risk") <- function() combined_risk

    return(fit)
}
