gamboostLSS_intern <- function(..., fun = c("check", "do_trace",
                                            "rescale_weights")) {

    fun <- match.arg(fun)
    do.call(fun, list(...))
}
