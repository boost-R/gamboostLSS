## helper functions

check <- function(what, what_char, names) {

    errormsg <- paste(sQuote(what_char), " can be either a scalar, a (named) vector or a (named) list",
                      " of ", what_char, " values with same names as ",  sQuote("families"), "in ",
                      sQuote("boost_control"))

    if (is.list(what)) {
        if (is.null(names(what)) && length(what) == length(names))
            names(what) <- names
        if (!all(names(what) %in% names) ||
            length(unique(names(what))) != length(names))
            stop(errormsg)
        what <- what[names] ## sort in order of families
        what <- unlist(what)
    } else {
        if(length(what) != 1 && length(what) != length(names))
            stop(errormsg)
        if (length(what) == 1) {
            what <- rep(what, length(names))
            names(what) <- names
        } else {
            if (is.null(names(what)))
                names(what) <- names
            if (!all(names(what) %in% names))
                stop(errormsg)
            what <- what[names] ## sort in order of families
        }
    }

    return(what)
}

## extract data from a gamboostLSS model (used in plot_PI)
get_data <- function(x, which = NULL) {
    data <- attr(x, "data")
    if (length(data) != 0) {
        data <- data[, which, drop = FALSE]
    } else {
        data <- try(sapply(which, get, env = parent.frame(2)),
                    silent = TRUE)
        if (inherits(data, "try-error"))
                stop("No data set found.")
        data <- as.data.frame(data)
    }
    return(data)
}

## extract family name from a gamboostLSS model (used in plot_PI)
get_families_name <- function(x) {
    attr(attr(x, "families"), "name")
}

## extract family name from a gamboostLSS model (used in plot_PI)
get_qfun <- function(x) {
    qfun <- attr(attr(x, "families"), "qfun")
    if (is.null(qfun))
        stop("Currently not implemented for this family")
    return(qfun)
}

## return mean or (first) modus of a vector depending on its class
mean_mod <- function(x) {
    if (is.numeric(x))
        return(mean(x, na.rm = TRUE))
    ## else compute and return modus
    if (is.character(x) || is.factor(x))
        return(names(which.max(table(x)))[1])
    stop("not implemented for data type ", class(x))
}

## helper function copied from mboost_2.2-3
rescale_weights <- function(w) {
    if (max(abs(w - floor(w))) < sqrt(.Machine$double.eps))
        return(w)
    return(w / sum(w) * sum(w > 0))
}

## helper function in a modified version based on mboost_2.2-3
## print trace of boosting iterations
do_trace <- function(current, mstart, risk,
                     linebreak = options("width")$width / 2, mstop = 1000) {
    current <- current - mstart
    if (current != mstop) {
        if ((current - 1) %/% linebreak == (current - 1) / linebreak) {
            mchr <- formatC(current + mstart, format = "d",
                            width = nchar(mstop) + 1, big.mark = "'")
            cat(paste("[", mchr, "] ",sep = ""))
        } else {
            if ((current %/% linebreak != current / linebreak)) {
                cat(".")
            } else {
                cat(" -- risk:", risk[current + mstart], "\n")
            }
        }
    } else {
        cat("\nFinal risk:", risk[current + mstart], "\n")
    }
}
