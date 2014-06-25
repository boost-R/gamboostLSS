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

    pred_vars <- lapply(x, extract, what = "variable.names")
    pred_vars <- unique(unlist(pred_vars))
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
