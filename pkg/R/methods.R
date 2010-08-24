### Methods

"[.mboostLSS" <- function(x, i, return = TRUE, ...) {
    stopifnot((length(i) == 1 | length(i) == length(x)) && i > 0)
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


risk <- function(object, ...)
    UseMethod("risk", object)

risk.mboostLSS <- function(object, merge = FALSE, parameter = names(object), ...){
    if (is.character(parameter))
        parameter <- extract_parameter(object, parameter)

    lo <- length(unique(mstop(object)))
    if (merge && lo == 1){
        RES <- sapply(parameter, function(i, object)  object[[i]]$risk(),
                      object = object)
        RES <- as.vector(t(RES))
        class(RES) <- object[[1]]$control$risk
        return(RES)
    }
    if (merge && lo != 1)
        warning(sQuote("merge = TRUE"), " not possible for multiple mstops")
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

selected.mboostLSS <- function(object, parameter = names(object), ...){
    if (is.character(parameter))
        parameter <- extract_parameter(object, parameter)
    RET <- lapply(parameter, function(i, object)  selected(object[[i]]),
                  object = object)
    names(RET) <- names(object)[parameter]
    if (length(RET) == 1)
        RET <- RET[[1]]
    return(RET)
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


fitted.mboostLSS <- function(object, parameter = names(object), ...){
    if (is.character(parameter))
        parameter <- extract_parameter(object, parameter)
    sapply(parameter, function(i, mod, ...) fitted(mod[[i]], ...),
           mod = object, ...)
}


predict.mboostLSS <- function(object, newdata = NULL,
                              type = c("link", "response", "class"),
                              which = NULL,
                              aggregate = c("sum", "cumsum", "none"),
                              parameter = names(object), ...) {
    if (is.character(parameter))
        parameter <- extract_parameter(object, parameter)
    sapply(parameter, function(i, mod, ...)
             predict(mod[[i]], newdata = newdata, type = type, which = which,
                     aggregate = aggregate, ...),
           mod = object, ...)
}







### helpers
do_trace <- function(current, mstart, mstop, risk,
                     linebreak = options("width")$width/2)
    mboost:::do_trace(m = current, mstop = mstart, risk = risk,
                      step = linebreak, width = mstop)


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
