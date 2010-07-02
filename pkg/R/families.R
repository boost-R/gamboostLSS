###
# Negative Binomial LSS Family

NBinomialMu <- function(off = list(sigma = 1)) {
    sigma <- off[["sigma"]]
    loss <- function(sigma, y, f)
        -(lgamma(y + sigma) - lgamma(sigma) - lgamma(y + 1) + sigma * log(sigma)
          - sigma * log(exp(f) + sigma) + y * f - y * log(exp(f) + sigma))
    risk <- function(y, f, w = 1){
        RET <- sum(w * loss(y = y, f = f, sigma = sigma))
        return(RET)
    }
    ngradient <- function(y, f, w = 1) {
        RET <- y - (y + sigma)/(exp(f) + sigma) * exp(f)
        return(RET)
    }
    Family(ngradient = ngradient, risk = risk, loss = loss,
           response = function(f) exp(f),
           check_y = function(y) {
               stopifnot(all.equal(unique(y - floor(y)), 0))
               y
           }, name = "Negative Negative-Binomial Likelihood (mu)")
}

NBinomialSigma <- function(off = list(mu = 0)) {
    # this family boosts sigma therefore f is sigma
    mu <- off[["mu"]]
    loss <- function(mu, y, f)
        -(lgamma(y + exp(f)) - lgamma(exp(f)) - lgamma(y + 1) + exp(f) * f -
          exp(f) * log(mu + exp(f)) + y * log(mu) - y * log(mu + exp(f)))

    offset <- function(y,w) 0
    risk <- function(y, f, w = 1){
        RET <- sum(w * loss(y = y, f = f, mu= mu))
        return(RET)
    }
    ngradient <- function(y, f, w = 1) {       # f is sigma !
        RET <- exp(f)*(digamma(y +exp(f)) - digamma(exp(f)) + log(exp(f)) + 1 -
                       log(mu +exp(f)) - (exp(f) + y)/(mu +exp(f)))
        return(RET)
    }
    Family(ngradient = ngradient,offset=offset, risk = risk, loss = loss,
           response = function(f) exp(f),
           check_y = function(y) {
               stopifnot(all.equal(unique(y - floor(y)), 0))
               y
           }, name = "Negative Negative-Binomial Likelihood (sigma)")
}


NBinomialLSS <- function(mu = 0, sigma = 1){
    RET <- list(mu = NBinomialMu(off = list(sigma = sigma)),
                sigma = NBinomialSigma(off = list(mu = mu)))
    class(RET) <- "families"
    return(RET)
}



###
# T-distribution LSS Family

StudentTMu <- function(off = list(df=1, sigma=1)) {
    df <- off[["df"]]
    sigma <- off[["sigma"]]

    loss <- function(df, sigma,y, f){
        -1 * (lgamma((df+1)/2) - log(sigma) - lgamma(1/2) - lgamma(df/2) - 0.5 *
              log(df) - (df+1)/2 * log(1 + (y-f)^2 / (df * sigma^2)))
    }

    offset <- function(y,w) 0
    risk <- function(y, f, w = 1){
        sum(w * loss(y = y, f = f, df = df, sigma = sigma))
    }
    ngradient <- function(y, f, w = 1) {
        (df+1)*(y-f)/(df*sigma^2 +(y-f)^2)
     }

    Family(ngradient = ngradient, risk = risk, loss = loss,
           offset=offset,
           response = function(f) f,
           check_y = function(y) {
               if (!is.numeric(y) || !is.null(dim(y)))
                   stop("response is not a numeric vector ")
               y
           },  name = "Student's t-distribution: mu (id link)")
}

StudentTSigma <- function(off = list(df=1, mu=0)) {
    mu <- off[["mu"]]
    df <- off[["df"]]

    loss <- function(df, mu, y, f){
        -1 * (lgamma((df+1)/2) - f - lgamma(1/2) - lgamma(df/2) - 0.5 * log(df) -
              (df+1)/2 * log(1 + (y-mu)^2 / (df * exp(2 * f))))
    }

    offset <- function(y,w) 1
    risk <- function(y, f, w = 1){
        sum(w * loss(y = y, f = f, df = df, mu = mu))
    }
    ngradient <- function(y, f, w = 1) {
        (-1 + ((df+1)*(y-mu)^2)/(df*exp(2*f) + (y-mu)^2/exp(2*f)))
    }

    Family(ngradient = ngradient, risk = risk, loss = loss,
           offset=offset,
           response = function(f) exp(f),
           check_y = function(y) {
               if (!is.numeric(y) || !is.null(dim(y)))
                   stop("response is not a numeric vector ")
               y
           },  name = "Student's t-distribution: sigma (log link)")
}

StudentTDf <- function(off = list(sigma=1, mu=0)) {
    mu <- off[["mu"]]
    sigma <- off[["sigma"]]

    loss <- function(sigma, mu,y, f){
        -1 * (lgamma((exp(f)+1)/2) - log(sigma) - lgamma(1/2) - lgamma(exp(f)/2) -
              0.5*f - (exp(f)+1)/2 * log(1 + (y-mu)^2 / (exp(f)*sigma^2)))
    }

    offset <- function(y,w) 1
    risk <- function(y, f, w = 1){
        sum(w * loss(y = y, f = f, sigma = sigma, mu = mu))
    }
    ngradient <- function(y, f, w = 1) {
        exp(f)/2 * (digamma((exp(f)+1)/2)-digamma(exp(f)/2)) - 0.5 -
            (exp(f)/2 * log(1+ (y-mu)^2/(exp(f)*sigma^2)) -
             (y-mu)^2/(sigma^2 + (y-mu)^2/exp(f)) * (exp(-f) +1)/2 )
    }

    Family(ngradient = ngradient, risk = risk, loss = loss,
           offset=offset,
           response = function(f) exp(f),
           check_y = function(y) {
               if (!is.numeric(y) || !is.null(dim(y)))
                   stop("response is not a numeric vector ")
               y
           },  name = "Student's t-distribution: df (log link)")
}



StudentTLSS <- function(mu = 0, sigma=1, df=1){
    RET <- list(mu = StudentTMu(off = list(sigma = sigma, df = df)),
                sigma = StudentTSigma(off = list(mu = mu, df = df)),
                df = StudentTDf(off = list(mu = mu, sigma = sigma)))
    class(RET) <- "families"
    return(RET)
}
