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

