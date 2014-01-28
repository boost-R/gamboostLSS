## Methods

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
    if (merge) {
        get_rsk <- function(i, object) {
            mmo <- max(mstop(object))
            rsk <- object[[i]]$risk()
            if (length(rsk) != mmo) {
                rsk <- c(rsk, rep(NA, mmo - length(rsk)))
            }
            rsk
        }
        RES <- sapply(parameter, get_rsk,
                      object = object)
        RES <- as.vector(t(RES))
        names(RES) <- rep(names(parameter), mstop(object)[1])
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
    RET <- lapply(parameter, function(i, object)
                               mboost::selected(object[[i]]),
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

plot_PI <- function(x, which, pi = 0.9, newdata = NULL,
                    main = "Marginal Prediction Intervals",
                    xlab = NULL, ylab = NULL, ...) {

    ## use check for families: currently only working for GaussianLSS
    avail <- c("Gaussian")
    if (is.null(get_families_name(x)) ||
        !(get_families_name(x) %in% avail))
        stop("currently only implemented for the following families:\n",
             paste(avail, collapse = ", "))

    if (length(which) != 1 || !is.character(which))
        stop("Please specify the variable name of the variable to be plotted")

    pred_vars <- lapply(x, extract, what = "variable.names")
    pred_vars <- unique(unlist(pred_vars))

    if (is.null(newdata)) {
        tmp <- get_data(x)[, pred_vars]
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

    if (is.null(xlab))
        xlab <- which
    if (is.null(ylab))
        ylab <- "prediction"

    plot(get_data(x)[, which], x[[1]]$response, pch = 20,
         col = rgb(0.5, 0.5, 0.5, 0.5),
         xlab = xlab, ylab = ylab, main = main)

    ## was ist denn bei anderen Verteilungen als der GaussianLSS?
    fm1 <- as.formula(paste(names(x)[1], " ~ ", which))
    fm2 <- as.formula(paste(names(x)[1], " + qnorm(1 - (1 - pi)/2) * ",
                            names(x)[2], " ~ ", which))
    fm3 <- as.formula(paste(names(x)[1], " + qnorm((1 - pi)/2) * ",
                            names(x)[2], " ~ ", which))

    lines(fm1, data = newdata)
    lines(fm2, data = newdata, lty = "dotted")
    lines(fm3, data = newdata, lty = "dotted")
}

get_data <- function(x) {
    attr(x, "data")
}

get_families_name <- function(x) {
    attr(attr(x, "families"), "name")
}

mean_mod <- function(x) {
    if (is.numeric(x))
        return(mean(x, na.rm = TRUE))
    ## else compute and return modus
    if (is.character(x) || is.factor(x))
        return(names(which.max(table(x))))
    stop("not implemented for data type ", class(x))
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

update.mboostLSS <- function(object, weights, oobweights = NULL,
                             risk = NULL, mstop = NULL, ...) {
    attr(object, "update")(weights = weights, oobweights = oobweights,
                           risk = risk, mstop = mstop, ...)
}

## generic version of model.weights (see stats::model.weights)
model.weights <- function(x, ...)
    UseMethod("model.weights")

model.weights.default <- function(x, ...)
    stats::model.weights(x)

model.weights.mboostLSS <- function(x, ...)
    attr(x, "(weights)")


################################################################################
### helpers

## print trace
do_trace <- function(current, mstart, mstop, risk,
                     linebreak = options("width")$width/2)
    mboost:::do_trace(m = current, mstop = mstart, risk = risk,
                      step = linebreak, width = mstop)


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
