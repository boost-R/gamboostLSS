## Methods

"[.mboostLSS" <- function(x, i, return = TRUE, ...) {
    stopifnot(length(i) == 1 | length(i) == length(x))
    attr(x, "subset")(i)
    if (return) return(x)
    invisible(NULL)
}

coef.mboostLSS <- function(object, which = NULL,
                           aggregate = c("sum", "cumsum", "none"),
                           parameter = names(object), ...){
    if (is.character(parameter))
        parameter <- extract_parameter(object, parameter)
    RET <- lapply(parameter, function(i, object)
        coef(object[[i]], which = which,
             aggregate = aggregate, ...),
        object = object)
    if (length(RET) == 1)
        RET <- RET[[1]]
    return(RET)
}

coef.glmboostLSS <- function(object, which = NULL,
                             aggregate = c("sum", "cumsum", "none"),
                             off2int = FALSE, parameter = names(object), ...){
    coef.mboostLSS(object, which = which, aggregate = aggregate,
                   parameter = parameter, off2int = off2int, ...)
}

risk.mboostLSS <- function(object, merge = FALSE, parameter = names(object), ...){
    if (is.character(parameter))
        parameter <- extract_parameter(object, parameter)
    
    lo <- length(unique(mstop(object)))
    if (merge) {
        get_rsk <- function(i, object) {
          mmo <- max(mstop(object)) + 1
          rsk <- object[[i]]$risk()
            if (length(rsk) != mmo) {
              if (length(rsk) > mmo)
                stop("risk cannot contain more entries than mstop")
                rsk <- c(rsk, rep(NA, mmo - length(rsk)))
            }
            rsk
        }
        RES <- sapply(parameter, get_rsk,
                      object = object)
        RES <- as.vector(t(RES))
        names(RES) <- rep(names(parameter), mstop(object)[1] + 1)
        ## drop unwanted NAs
        if (lo != 1)
            RES <- RES[!is.na(RES)]
        class(RES) <- object[[1]]$control$risk
        return(RES)
    }
    RES <- lapply(parameter, function(i, object)  object[[i]]$risk(),
                  object = object)
    class(RES) <- object[[1]]$control$risk
    return(RES)
}

mstop.mboostLSS <- function(object, parameter = names(object), ...){
    if (is.character(parameter))
        parameter <- extract_parameter(object, parameter)
    RET <- sapply(parameter, function(i, object)  object[[i]]$mstop(),
                  object = object)
    names(RET) <- names(object)[parameter]
    if (length(RET) == 1)
        RET <- RET[[1]]
    # change mstop for noncyclical fitting to scalar value with attribute of partitions
    if (inherits(object, "nc_mboostLSS")) {
        partitions <- RET
        RET <- sum(RET)
        attr(RET, "partitions") <- partitions
    }
    
    return(RET)
}

mstop.oobag <- function(object, parameter = names(object), ...){
    if (is.character(parameter))
        parameter <- extract_parameter(object, parameter)
    RET <- sapply(parameter, function(i, object)  which.min(object[[i]]),
                  object = object)
    names(RET) <- names(object)[parameter]
    if (length(RET) == 1)
        RET <- RET[[1]]
    return(RET)
}

selected.mboostLSS <- function(object, merge = FALSE, parameter = names(object),
                               ...){
    if (is.character(parameter))
        parameter <- extract_parameter(object, parameter)
    
    #merge is different for noncyclical fitting
    if (merge) {
        if (inherits(object, "nc_mboostLSS")){
            
            #get the names of parameter selected in each iteration (drop initial offset risk values)
            RET <- names(attr(object, "combined_risk")())[-seq_along(parameter)]
            names(RET) <- RET #set the names of the vector as we will overwrite the values.
            
            #overwrite names in the vector with the selected BLs in the correct order
            for(p in names(parameter)){
                RET[RET == p] <- object[[p]]$xselect()
            }
            mode(RET) = "numeric" #ensure numeric values -> as.numeric drops the names
         
            return(RET)
        }
        else {
            get_sel <- function(i, object) {
                mmo <- max(mstop(object))
                sel <- object[[i]]$xselect()
                if (length(sel) != mmo) {
                    sel <- c(sel, rep(NA, mmo - length(sel)))
                }
                sel
            }
            
            RET <- sapply(parameter, get_sel,
                          object = object)
            
            RET <- as.vector(t(RET))
            names(RET) <- rep(names(parameter), mstop(object)[1])
            lo <- length(unique(mstop(object)))
            if (lo != 1)
                RET <- RET[!is.na(RET)]
            return(RET)
        }
        
    }
    else {
        RET <- lapply(parameter, function(i, object)
            selected(object[[i]]),
            object = object)
        names(RET) <- names(object)[parameter]
        if (length(RET) == 1)
            RET <- RET[[1]]
        return(RET)
    }
    
}

selected.stabsel_mboostLSS <- function(object, parameter = NULL, ...) {
    ## extract parameters
    param_count <- table(gsub(".*\\.", "", rownames(object$phat)))
    ## set up a named list for the results
    res <- vector("list", length(param_count))
    names(res) <- names(param_count)
    ## subtract the number of base-learners of the first components
    param_count <- c(0, cumsum(param_count))
    for (i in 1:length(res)) {
        idx <- gsub(".*\\.", "", names(object$selected)) == names(res)[[i]]
        res[[i]] <- object$selected[idx] - param_count[i]
    }
    if (!is.null(parameter))
        res <- res[parameter]
    if (length(res) == 1)
        res <- res[[1]]
    return(res)
}


plot.glmboostLSS <- function(x, main = names(x), parameter = names(x),
                             off2int = FALSE, ...){
    if (is.character(parameter))
        parameter <- extract_parameter(x, parameter)
    lapply(parameter, function(i, x, main, off2int,  ...)
        plot(x[[i]], main = main[[i]], off2int = off2int,  ...),
        x = x, main = main, off2int = off2int, ...)
    invisible(coef(x, aggregate = "cumsum", off2int = off2int))
}


plot.gamboostLSS <- function(x, main = names(x), parameter = names(x), ...){
    if (is.character(parameter))
        parameter <- extract_parameter(x, parameter)
    RET <- lapply(parameter, function(i, x, main, ...)
        plot(x[[i]], main = main[[i]], ...),
        x = x, main = main, ...)
    if (any(sapply(RET, class) == "trellis")) {
        return(RET)
    } else {
        invisible(RET)
    }
}

plot.predint <- function(x, main = "Marginal Prediction Interval(s)",
                         xlab = NULL, ylab = NULL, lty = c("solid", "dashed"),
                         lcol = c("black", "black"), log = "", ...) {
    
    pi <- attr(x, "pi")
    which <- attr(x, "which")
    rawdata <- attr(x, "rawdata")
    
    if (length(lty) != length(pi) + 1)
        lty <- c(lty, rep(tail(lty, 1), (length(pi) + 1) - length(lty)))
    if (length(lcol) != length(pi) + 1)
        lcol <- c(lcol, rep(tail(lcol, 1), (length(pi) + 1) - length(lcol)))
    
    if (is.null(xlab))
        xlab <- which
    if (is.null(ylab))
        ylab <- "prediction"
    
    plot(rawdata$x, rawdata$y, pch = 20,
         col = rgb(0.5, 0.5, 0.5, 0.5),
         xlab = xlab, ylab = ylab, main = main,
         log = log, ...)
    
    lines(x[, which], x$"Prediction (Median)",
          lty = lty[1], col = lcol[1], ...)
    
    for (i in seq_along(pi)) {
        lines(x[, which], x[, paste0(pi[i] * 100, "% PI (lower)")],
              lty = lty[i + 1], col = lcol[i + 1], ...)
        lines(x[, which], x[, paste0(pi[i] * 100, "% PI (upper)")],
              lty = lty[i + 1], col = lcol[i + 1], ...)
    }
}


PI <- predint <- function(x, which, pi = 0.9, newdata = NULL, ...) {
    qfun <- get_qfun(x)
    
    if (length(which) != 1 || !is.character(which))
        stop("Please specify the variable for the marginal prediction interval.")
    
    var <- get_data(x, which = which)
    if (ncol(var) > 1 || is.factor(var))
        stop("Prediction intervals are currently only implemented for ",
             "base-learners of one numeric variable")
    
    pred_vars <- lapply(x, extract, what = "variable.names")
    pred_vars <- unique(trimws(unlist(strsplit(unlist(pred_vars), ",")))) # fixing bug #55
    if ("(Intercept)" %in% pred_vars)
        pred_vars <- pred_vars[pred_vars != "(Intercept)"]
    
    if (is.null(newdata)) {
        tmp <- get_data(x, which = pred_vars)
        i <- grepl(which, names(tmp))
        if (sum(i) != 1)
            stop(sQuote("which"), " is misspecified")
        newdata <- data.frame(x1 = seq(min(tmp[, i]),
                                       max(tmp[, i]), length = 150))
        colnames(newdata) <- names(tmp)[i]
        newdata[, names(tmp)[!i]] <- lapply(tmp[, !i], mean_mod)
    } else {
        i <- grepl(which, names(newdata))
        if (sum(i) != 1)
            stop(sQuote("which"), "is misspecified")
        ## check if data is ok, else give a warning
        if (nrow(unique(newdata[, !i])) != 1)
            warning("All variables but", sQuote("which"), "should be constant")
    }
    
    newdata[, names(x)] <- predict(x, newdata = newdata, type = "response")
    newdata$"Prediction (Median)" <- do.call(qfun, args = c(p = 0.5,
                                                            newdata[, names(x)]))
    
    for (i in seq_along(pi)) {
        newdata[, paste0(pi[i] * 100, "% PI (lower)")] <- do.call(qfun,
                                                                  args = c(p = (1 - pi[i])/2, newdata[, names(x)]))
        newdata[, paste0(pi[i] * 100, "% PI (upper)")] <- do.call(qfun,
                                                                  args = c(p = 1 - (1 - pi[i])/2, newdata[, names(x)]))
    }
    # drop predictions of parameters
    newdata <- newdata[, !names(newdata) %in% names(x)]
    
    class(newdata) <- c("predint", "data.frame")
    attr(newdata, "pi") <- pi
    attr(newdata, "which") <- which
    attr(newdata, "rawdata") <- data.frame(x = get_data(x, which = pred_vars)[, which],
                                           y = x[[1]]$response)
    return(newdata)
}


print.mboostLSS <- function(x, ...){
    cl <- match.call()
    cat("\n")
    cat("\t LSS Models fitted via Model-based Boosting\n")
    cat("\n")
    if (!is.null(attr(x, "call")))
        cat("Call:\n", deparse(attr(x, "call")), "\n\n", sep = "")
    m <- mstop(x)
    if (inherits(x, "nc_mboostLSS"))
        m <- attr(m, "partitions")
    cat("Number of boosting iterations (mstop): ",
        paste(names(x), m, sep = " = ", collapse = ", "), "\n")
    nus <- sapply(x, function(xi) xi$control$nu)
    cat("Step size: ",
        paste(names(nus), nus, sep = " = ", collapse = ", "), "\n\n")
    cat("Families:\n")
    lapply(x, function(xi) show(xi$family))
    invisible(x)
}


fitted.mboostLSS <- function(object, parameter = names(object), ...){
    if (is.character(parameter))
        parameter <- extract_parameter(object, parameter)
    myApply(parameter, function(i, mod, ...) fitted(mod[[i]], ...),
            mod = object, ...)
}


predict.mboostLSS <- function(object, newdata = NULL,
                              type = c("link", "response", "class"),
                              which = NULL,
                              aggregate = c("sum", "cumsum", "none"),
                              parameter = names(object), ...) {
    if (is.character(parameter))
        parameter <- extract_parameter(object, parameter)
    myApply(parameter, function(i, mod, ...)
        predict(mod[[i]], newdata = newdata, type = type, which = which,
                aggregate = aggregate, ...),
        mod = object, ...)
}

update.mboostLSS <- function(object, weights, oobweights = NULL,
                             risk = NULL, trace = NULL, mstop = NULL, ...) {
    attr(object, "update")(weights = weights, oobweights = oobweights,
                           risk = risk, trace = trace, mstop = mstop, ...)
}

## generic version of model.weights (see stats::model.weights)
model.weights <- function(x, ...)
    UseMethod("model.weights")

model.weights.default <- function(x, ...)
    stats::model.weights(x)

model.weights.mboostLSS <- function(x, ...)
    attr(x, "(weights)")

### summary function based on print.mboostLSS() and summary.mboost()
summary.mboostLSS <- function(object, ...) {
    cat("\n")
    cat("\t LSS Models fitted via Model-based Boosting\n")
    cat("\n")
    if (!is.null(attr(object, "call")))
        cat("Call:\n", deparse(attr(object, "call")), "\n\n", sep = "")
    m <- mstop(object)
    if (inherits(object, "nc_mboostLSS"))
        m <- attr(m, "partitions")
    cat("Number of boosting iterations (mstop): ",
        paste(names(object), m, sep = " = ", collapse = ", "), "\n")
    nus <- sapply(object, function(xi) xi$control$nu)
    cat("Step size: ",
        paste(names(nus), nus, sep = " = ", collapse = ", "), "\n\n")
    
    cat("Families:\n")
    lapply(object, function(xi) show(xi$family))
    
    if (inherits(object, "glmboostLSS")) {
        cat("Coefficients:\n")
        cf <- coef(object, off2int = TRUE)
        for (i in 1:length(cf)) {
            cat("Parameter ", names(cf)[i], ":\n", sep = "")
            print(cf[[i]])
            cat("\n")
        }
    }
    
    if (!all(is_null <- sapply(selected(object), is.null))) {
        cat("Selection frequencies:\n")
        for (i in 1:length(object)) {
            cat("Parameter ", names(object)[i], ":\n", sep = "")
            if (is_null[i]){
                print(NULL)
                next
            }
            nm <- variable.names(object[[i]])
            selprob <- tabulate(selected(object[[i]]), nbins = length(nm)) /
                length(selected(object[[i]]))
            names(selprob) <- names(nm)
            selprob <- sort(selprob, decreasing = TRUE)
            selprob <- selprob[selprob > 0]
            print(selprob)
        }
    }
    invisible(object)
}

stabsel.mboostLSS <- function(x, cutoff, q, PFER,
                              mstop = NULL,
                              folds = subsample(model.weights(x), B = B),
                              B = ifelse(sampling.type == "MB", 100, 50),
                              assumption = c("unimodal", "r-concave", "none"),
                              sampling.type = c("SS", "MB"),
                              papply = mclapply, verbose = TRUE, FWER, eval = TRUE, ...) {
    
    cll <- match.call()
    p <- sum(sapply(x, function(obj) length(variable.names(obj))))
    
    if(inherits(x, "FDboostLSS")) {
      if(is.null(x[[1]]$ydim)){
        n <- length(attr(x, "(weights)")) # scalar response
      }else{
        n <- x[[1]]$ydim[1] # functional reponse 
        # correct the wrong default folds if necessary
        if(nrow(folds) == length(model.weights(x))){
          folds <- subsample(rep(1, n), B = B)
        }
      }
    } else {
      n <- nrow(attr(x, "data"))
    }
    
    ## extract names of base-learners (and add paramter name)
    nms <- lapply(x, function(obj) variable.names(obj))
    nms <- lapply(1:length(nms), function(i)
        paste(nms[[i]], names(nms)[i], sep = "."))
    
    sampling.type <- match.arg(sampling.type)
    if (sampling.type == "MB")
        assumption <- "none"
    else
        assumption <- match.arg(assumption)
    
    
    if (inherits(x, "nc_mboostLSS")) {
        ## check mstop
        if (is.null(mstop))
            mstop <- sum(mstop(x))
        if (length(mstop) != 1 | mstop %% 1 != 0 | mstop < length(x)) {
            stop(sQuote("mstop"), " has to be an integer larger than ",
                 length(x))
        }
    }
    else {
        if (is.null(mstop))
            mstop <- mstop(x)
        mstop <- check(mstop, "mstop", names(x))
    }
    
    if (length(unique(mstop)) != 1)
        warning("Usually one should use the same ", sQuote("mstop"),
                " value for all components for stability selection.")
    
    if (verbose)
        cat("Run stabsel ")
    
    ## set mstop = 0 to speed things up
    x <- update(x, weights = model.weights(x), mstop = 0)
    
    ## define the fitting function (args.fitfun is not used but needed for
    ## compatibility with run_stabsel
    fit_model <- function(i, folds, q, args.fitfun) {
        if (verbose)
            cat(".")
        ## start by setting up model on subset and fit first q iterations
        mod <- update(x, weights = folds[, i], mstop = q)
        ## make sure dispatch works correctly
        class(mod) <- class(x)
        xs <- selected(mod)
        nsel <- length(mod)
        ## now update model until we obtain q different base-learners altogether
        for (m in (q+1):max(mstop)) {
            if (nsel >= q)
                break
            mstop(mod) <- m
            xs <- selected(mod)
            nsel <- sum(sapply(xs, function(selection)
                length(unique(selection))))
        }
        #this changes nothing for method = "cyclic" but fixes mstop for method = "noncyclic"
        mstop <- check(mstop, "mstop", names(x))
        ## complete paths
        if (any(sapply(xs, length) < mstop)) {
            for (j in 1:length(xs)) {
                
## <FIXME> What happens if component j was never selected, i.e. xs[[j]] = NULL?
## Can we use NA as proposed? We need to see what happens later.
                if (is.null(xs[[j]]))
                    xs[[j]][1] <- NA
                start <- length(xs[[j]]) + 1
                xs[[j]][start:mstop[j]] <- xs[[j]][1]
## </FIXME>   
                
            }
        }
        
        selected <- lapply(xs, unique)
        ret <- lapply(1:length(selected), function(i) {
            res <- logical(length(nms[[i]]))
            names(res) <- nms[[i]]
            res[selected[[i]]] <- TRUE
            res
        })
        ret <- unlist(ret)
        
        ## compute selection paths
        #merging for method cyclic
        if(!inherits(x, "nc_mboostLSS")){
            sequences <- lapply(1:length(xs), function(i) {
                res <- matrix(FALSE, nrow = length(nms[[i]]), ncol = mstop[[i]])
                rownames(res) <- nms[[i]]
                for (j in 1:mstop[[i]])
                    res[xs[[i]][j], j:mstop[[i]]] <- TRUE
                res
            })

## <FIXME> What is this error message about? No user will know what you mean. Please fix coding issue and remove stop() or 
## provide a relevant error message.            
            if (any(mstop < max(mstop)))
                stop("simply add the last column to the smaller matrices")
## </FIXME>
            
            ## now merge sequences
            for (i in 1:ncol(sequences[[1]])) {
                for (j in 1:length(sequences)) {
                    if (i == 1) {
                        if (j == 1) {
                            other_params <- rep(FALSE, sum(sapply(sequences, nrow)[-1]))
                            sequence <- matrix(c(sequences[[i]][, j],
                                                 other_params))
                        } else {
                            tmp <- unlist(lapply(sequences[1:j], function(x) x[, i]))
                            other_params <- rep(FALSE, sum(sapply(sequences,
                                                                  nrow)[-(1:j)]))
                            tmp <- c(tmp, other_params)
                            sequence <- cbind(sequence, tmp)
                        }
                    } else {
                        if (j < length(sequences)) {
                            tmp <- unlist(c(lapply(sequences[1:j], function(x) x[, i]),
                                            lapply(sequences[(j+1):length(sequences)],
                                                   function(x) x[, i - 1])))
                        } else {
                            tmp <- unlist(lapply(sequences[1:j], function(x) x[, i]))
                        }
                        sequence <- cbind(sequence, tmp)
                    }
                }
            }
        }
        else {
            sequence <- matrix(FALSE, nrow = p, ncol = mstop[1])
            rownames(sequence) <- unlist(nms)
            
            sel <- selected(mod, merge = TRUE)
            
            for (i in names(mod)) {
                varnames <- variable.names(mod[[i]])
                for(j in seq_along(varnames)){
                    pos <- which(names(sel) == i & sel == j)
                    if(length(pos) > 0)
                        sequence[paste(varnames[j], i, sep = "."), min(pos):mstop[1]] <- TRUE
                }
                
            }
            
        }
        
        colnames(sequence) <- 1:ncol(sequence)
        ret <- list(selected = ret, path = sequence)
        ## was mstop to small?
        attr(ret, "violations") <- ifelse(sum(ret$selected) < q, TRUE, FALSE)
        return(ret)
    }
    
    ret <- run_stabsel(fitter = fit_model, args.fitter = list(),
                       n = n, p = p, cutoff = cutoff, q = q,
                       PFER = PFER, folds = folds, B = B, assumption = assumption,
                       sampling.type = sampling.type, papply = papply,
                       verbose = verbose, FWER = FWER, eval = eval,
                       names = unlist(nms), ...)
    
    if (verbose)
        cat("\n")
    
    if (!eval)
        return(ret)
    
    violations <- FALSE
    if (!is.null(attr(ret, "violations")))
        violations <- attr(ret, "violations")
    
    if (any(violations))
        warning(sQuote("mstop"), " too small in ",
                sum(violations), " of the ", length(violations),
                " subsampling replicates to select ", sQuote("q"),
                " base-learners; Increase ", sQuote("mstop"),
                " bevor applying ", sQuote("stabsel"))
    
    ret$call <- cll
    ret$call[[1]] <- as.name("stabsel")
    class(ret) <- c("stabsel_mboostLSS", "stabsel")
    ret
}

################################################################################
### helpers

## extract parameter index from mboostLSS object x
extract_parameter <- function(x, parameter) {
    idx <- sapply(parameter, function(w) {
        wi <- grep(w, names(x), fixed = TRUE)
        if (length(wi) > 0) return(wi)
        return(NA)
    })
    if (any(is.na(idx)))
        warning(paste(parameter[is.na(idx)], collapse = ","), " not found")
    parameter <- idx
}

## function for weighted sd
weighted.sd <- function(x, w, ...) {
    if (missing(w))
        w <- rep(1, length(x))
    m <- weighted.mean(x, w, ...)
    var <- weighted.mean((x - m)^2, w, ...) * sum(w) / (sum(w) - 1)
    return(sqrt(var))
}

## weighted median
weighted.median <- function (x, w = 1, na.rm = FALSE) {
    if (length(w) == 1)
        w <- rep(w, length(x))
    
    ## remove observations with zero weights
    x <- x[w != 0]
    w <- w[w != 0]
    
    ## remove NAs if na.rm = TRUE
    if (na.rm) {
        keep <- !is.na(x) & !is.na(w)
        x <- x[keep]
        w <- w[keep]
    } else {
        if (any(is.na(x)) | any(is.na(w)))
            return(NA)
    }
    
    ## sort data and weights
    ind <- order(x)
    x <- x[ind]
    w <- w[ind]
    
    ## first time that fraction of weights is above 0.5
    ind1 <- min(which(cumsum(w)/sum(w) > 0.5))
    
    ## first time that fraction of weights is below 0.5
    ind2 <- ifelse(ind1 == 1, 1, max(which(cumsum(w)/sum(w) <= 0.5)))
    
    ## if sum of weights is an even integer
    if(sum(w) %% 1 == 0 && sum(w) %% 2 == 0)
        return(mean(c(x[ind1], x[ind2])))
    
    ## else return
    return(max(c(x[ind1], x[ind2])))
}
