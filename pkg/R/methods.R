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


plot.glmboostLSS <- function(x, ...){
    lapply(x, plot, ...)
}


print.mboostLSS <- function(x, ...){
    lapply(x, print, ...)
}
