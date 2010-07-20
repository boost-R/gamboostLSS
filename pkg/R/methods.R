### Methods

risk <- function(object, ...)
    UseMethod("risk", object)

risk.mboostLSS <- function(object, merge = FALSE, ...){
    if (merge){
        RES <- sapply(object, function(x) x$risk())
        RES <- as.vector(t(RES))
        class(RES) <- object[[1]]$control$risk
        return(RES)
    } else {
        RES <- lapply(object, function(x) x$risk())
        class(RES) <- object[[1]]$control$risk
        return(RES)
    }
}

coef.mboostLSS <- function(object, ...)
    lapply(object, coef, ...)

mstop.mboostLSS <- function(object, ...){
    RET <- sapply(object, function(x) x$mstop())[length(object)]
    names(RET) <- NULL
    RET
}

mstop.oobag <- function(object, ...){
    RET <- sapply(object, which.min)[length(object)]
    names(RET) <- NULL
    RET
}

"[.mboostLSS" <- function(x, i, return = TRUE, ...) {
    stopifnot(length(i) == 1 && i > 0)
    attr(x, "subset")(i)
    if (return) return(x)
    invisible(NULL)
}


selected.mboostLSS <- function(object, ...){
    lapply(object, selected)
}


plot.glmboostLSS <- function(x, main = names(x), parameter = names(x),
                             off2int = FALSE, ...){
    if (is.character(parameter)) {
        idx <- sapply(parameter, function(w) {
            wi <- grep(w, names(x), fixed = TRUE)
            if (length(wi) > 0) return(wi)
            return(NA)
        })
        if (any(is.na(idx)))
            warning(paste(parameter[is.na(idx)], collapse = ","), " not found")
        parameter <- idx
    }
    lapply(parameter, function(i, x, main, off2int,  ...)
                    plot(x[[i]], main = main[[i]], off2int = off2int,  ...),
           x = x, main = main, off2int = off2int, ...)
    invisible(coef(x, aggregate = "cumsum", off2int = off2int))
}


plot.gamboostLSS <- function(x, main = names(x), parameter = names(x), ...){
    if (is.character(parameter)) {
        idx <- sapply(parameter, function(w) {
            wi <- grep(w, names(x), fixed = TRUE)
            if (length(wi) > 0) return(wi)
            return(NA)
        })
        if (any(is.na(idx)))
            warning(paste(parameter[is.na(idx)], collapse = ","), " not found")
        parameter <- idx
    }
    lapply(parameter, function(i, x, main, ...)
                    plot(x[[i]], main = main[[i]], ...),
           x = x, main = main, ...)
    invisible(x)
}


print.mboostLSS <- function(x, ...){
    cl <- match.call()
    cat("\n")
    cat("\t LSS Models fitted via Model-based Boosting\n")
    cat("\n")
    if (!is.null(attr(x, "call")))
        cat("Call:\n", deparse(attr(x, "call")), "\n\n", sep = "")
    cat("Number of boosting iterations: mstop =", mstop(x), "\n")
    cat("Step size: ", x[[1]]$control$nu, "\n\n")
    cat("Families:\n")
    lapply(x, function(xi) show(xi$family))
    invisible(x)
}
