###################################################################
## cross-validation (bootstrap, k-fold cv etc.) of empirical risk
## for boosting algorithms for gamLSS models

make.grid <- function(max, length.out = 10, min = NULL, log = TRUE) {
    if (is.null(min))
        min <- rep(1, length(max))
    if (length(min) == 1)
        min <- rep(min, length(max))
    if (length(length.out) == 1)
        length.out <- rep(length.out, length(max))
    if (length(length.out) != length(max))
        stop(sQuote("length.out"),
             " must be either scalar or a vector of the same length as ",
             sQuote("max"))
    if (length(min) != length(max))
        stop(sQuote("min"),
             " must be either scalar or a vector of the same length as ",
             sQuote("max"))

    if (log == TRUE) {
        min <- log(min)
        max <- log(max)
    }

    if (any(sapply(1:length(max), function(i) min[i] >= max[i])))
        stop("All min values must be smaller than the respectiv max value.")

    if (length(max) == 1) {
        RET <- seq(from = min, to = max, length.out = length.out)
        if (log == TRUE)
            RET <- exp(RET)
        ## round to get integer values
        RET <- round(RET)
        if (any(duplicated(RET))) {
            warning("Duplicates produced; Only unique values are returned")
            RET <- unique(RET)
        }
        return(RET)
    }

    RET <- lapply(1:length(max), function(i)
                  seq(from = min[i], to = max[i], length.out = length.out[i]))
    if (log == TRUE)
        RET <- lapply(RET, exp)
    ## round to get integer values
    RET <- lapply(RET, round)
    if (any(sapply(RET, duplicated))) {
        warning("Duplicates produced; Only unique values are returned")
        RET <- lapply(RET, unique)
    }

    if (!is.null(names(max)))
        names(RET) <- names(max)
    return(RET)
}


###
# cvrisk, adapted version from mboost (2.2-2)
cvrisk.mboostLSS <- function(object, folds = cv(model.weights(object)),
                             grid = make.grid(mstop(object)),
                             papply = mclapply,
                             fun = NULL, ...) {

    weights <- model.weights(object)
    if (any(weights == 0))
        warning("Zero weights in ", sQuote("object"))
    if (is.null(folds)) {
        folds <- rmultinom(25, length(weights), weights/sum(weights))
    } else {
        stopifnot(is.matrix(folds) && nrow(folds) == length(weights))
    }
    if (length(object) != length(grid))
        stop(sQuote("grid"),
             " must be a list of the same length as parameters in",
             sQuote("object"))

    oobrisk <- matrix(0, nrow = ncol(folds), ncol = length(grid))
    if (!is.null(fun))
        stopifnot(is.function(fun))

    ### WHAT ABOUT:
    ## fam_name <- object$family@name
    call <- deparse(attr(object, "call"))

    if (is.null(fun)) {
        dummyfct <- function(weights, oobweights) {
            ## make model with new weights and minimal mstop
            mod <- update(object, weights = weights, oobweights = oobweights,
                          risk = "oobag", mstop = lapply(grid, min))
            ## now we need to increase mstop (stupid or clever)
            x_grid <- expand.grid(grid)
            risks <- vector("numeric", nrow(x_grid))
            for (i in 1:nrow(x_grid)) {
                mod[x_grid[i, ]]
                rsk <- risk(mod, merge = TRUE)
                risks[i] <- rsk[length(rsk)]
            }
            return(risks)
        }
    }
    else {
        stop("currently not implemented")
        dummyfct <- function(weights, oobweights) {
            mod <- update(object, weights = weights, oobweights = oobweights,
                          risk = "oobag", mstop = lapply(grid, min))
            mod[max(grid)]
            ### make sure dispatch works correctly
            class(mod) <- class(object)
            fun(mod)
        }
    }

    OOBweights <- matrix(rep(weights, ncol(folds)), ncol = ncol(folds))
    OOBweights[folds > 0] <- 0
    oobrisk <- papply(1:ncol(folds),
        function(i) dummyfct(weights = folds[, i],
                             oobweights = OOBweights[, i]), ...)
    ## get errors if mclapply is used
    if (any(idx <- sapply(oobrisk, is.character)))
        stop(sapply(oobrisk[idx], function(x) x))
    if (!is.null(fun))
        return(oobrisk)
    oobrisk <- t(as.data.frame(oobrisk))
    oobrisk <- oobrisk / colSums(OOBweights)
    colnames(oobrisk) <- apply(expand.grid(grid), 1,
                               function(x) paste(x, collapse = ","))
    rownames(oobrisk) <- 1:nrow(oobrisk)
    ## attr(oobrisk, "risk") <- fam_name
    ## sowas wie "Normal distribution: mu(id link)"
    attr(oobrisk, "call") <- call
    attr(oobrisk, "mstop") <- grid
    attr(oobrisk, "type") <- ifelse(!is.null(attr(folds, "type")),
        attr(folds, "type"), "user-defined")
    class(oobrisk) <- "cvriskLSS"
    oobrisk
}

print.cvriskLSS <- function(x, ...) {
    #cat("\n\t Cross-validated", attr(x, "risk"), "\n\t",
    cat("\n\t Cross-validated risk\n\t",
              attr(x, "call"), "\n\n")
    print(colMeans(x))
    cat("\n\t Optimal number of boosting iterations:", mstop(x), "\n")
    return(invisible(x))
}

plot.cvriskLSS <- function(x, ylab = attr(x, "risk"),
                           xlab = "Number of boosting iterations",
                           ylim = range(x), main = attr(x, "type"), ...) {

    warning("This function is preliminary only")
    cm <- colMeans(x)
    plot(1:ncol(x), cm, ylab = ylab, ylim = ylim,
         type = "n", lwd = 2, xlab = xlab,
         main = main, axes = FALSE, ...)
    out <- apply(x, 1, function(y) lines(1:ncol(x),y, col = "lightgrey"))
    rm(out)
    ms <- which.min(cm)
    lines(c(ms, ms), c(min(c(0, ylim[1] * ifelse(ylim[1] < 0, 2, 0.5))), cm[ms]),
          lty = 2)
    lines(1:ncol(x), cm, type = "l")
    axis(1, at = 1:ncol(x), labels = colnames(x))
    axis(2)
    box()
}

mstop.cvriskLSS <- function(object, parameter = NULL, ...) {
    res <- unlist(expand.grid(attr(object, "mstop"))[which.min(colSums(object)),])
    if (!is.null(parameter)) {
        if(is.character(parameter))
            parameter <- extract_parameter(object, parameter)
        res <- res[parameter]
    }
    return(res)
}
