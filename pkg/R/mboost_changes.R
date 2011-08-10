## functions from package mboost that need to be changed in order to make
## gamboostLSS work with mboost <= 2.0.11

fitted.mboost <- function(object, ...) {
    args <- list(...)
    if (length(args) == 0) {
        ret <- object$fitted()
        names(ret) <- object$rownames
    } else {
        ret <- predict(object, newdata=NULL, ...)
        if (NROW(ret) == length(ret))
            rownames(ret) <- object$rownames
    }
    ret
}

selected <- function(object, ...)
    UseMethod("selected", object)

plot.glmboost <- function(x, main = deparse(x$call), col = NULL,
						  off2int = FALSE, ...) {

    cp <- coef(x, aggregate = "cumsum", off2int = off2int)
    ncp <- names(cp)
    cp <- matrix(unlist(cp), nrow = length(cp), byrow = TRUE)
    cf <- cp[, ncol(cp)]
    if (is.null(col))
        col <- hcl(h = 40, l = 50, c= abs(cf) / max(abs(cf)) * 490)
    matplot(t(cp), type = "l", lty = 1, xlab = "Number of boosting iterations",
            ylab = "Coefficients", main = main, col = col, ...)
    abline(h = 0, lty = 1, col = "lightgray")
    axis(4, at = cf, labels = ncp, las = 1)
}
