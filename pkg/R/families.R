###
# Negative Binomial LSS Family

NBinomialMu <- function(mu = NULL, sigma = NULL) {
    offset <- function(y, w)
        log(mu)
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
    Family(ngradient = ngradient, offset = offset, risk = risk, loss = loss,
           response = function(f) exp(f),
           check_y = function(y) {
               stopifnot(all.equal(unique(y - floor(y)), 0))
               y
           }, name = "Negative Negative-Binomial Likelihood: mu (log link)")
}

NBinomialSigma <- function(mu = NULL, sigma = NULL) {
    # this family boosts sigma therefore f is sigma
    loss <- function(mu, y, f)
        -(lgamma(y + exp(f)) - lgamma(exp(f)) - lgamma(y + 1) + exp(f) * f -
          exp(f) * log(mu + exp(f)) + y * log(mu) - y * log(mu + exp(f)))

    offset <- function(y,w)
        log(sigma)
    risk <- function(y, f, w = 1){
        RET <- sum(w * loss(y = y, f = f, mu= mu))
        return(RET)
    }
    ngradient <- function(y, f, w = 1) {       # f is sigma !
        RET <- exp(f)*(digamma(y +exp(f)) - digamma(exp(f)) + log(exp(f)) + 1 -
                       log(mu +exp(f)) - (exp(f) + y)/(mu +exp(f)))
        return(RET)
    }
    Family(ngradient = ngradient, offset = offset, risk = risk, loss = loss,
           response = function(f) exp(f),
           check_y = function(y) {
               stopifnot(all.equal(unique(y - floor(y)), 0))
               y
           }, name = "Negative Negative-Binomial Likelihood: sigma (log link)")
}


NBinomialLSS <- function(mu = 1, sigma = 1){
    if (sigma <= 0)
        stop(sQuote("sigma"), " must be greater than zero")
    RET <- list(mu = NBinomialMu(mu = mu, sigma = sigma),
                sigma = NBinomialSigma(mu = mu, sigma = sigma))
    class(RET) <- "families"
    return(RET)
}



###
# T-distribution LSS Family

StudentTMu <- function(mu = NULL, sigma = NULL, df = NULL) {
    loss <- function(df, sigma,y, f){
        -1 * (lgamma((df+1)/2) - log(sigma) - lgamma(1/2) - lgamma(df/2) - 0.5 *
              log(df) - (df+1)/2 * log(1 + (y-f)^2 / (df * sigma^2)))
    }

    offset <- function(y,w)
        mu
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

StudentTSigma <- function(mu = NULL, sigma = NULL, df = NULL) {
    loss <- function(df, mu, y, f){
        -1 * (lgamma((df+1)/2) - f - lgamma(1/2) - lgamma(df/2) - 0.5 * log(df) -
              (df+1)/2 * log(1 + (y-mu)^2 / (df * exp(2 * f))))
    }

    offset <- function(y,w)
        log(sigma)
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

StudentTDf <- function(mu = NULL, sigma = NULL, df = NULL) {
    loss <- function(sigma, mu,y, f){
        -1 * (lgamma((exp(f)+1)/2) - log(sigma) - lgamma(1/2) - lgamma(exp(f)/2) -
              0.5*f - (exp(f)+1)/2 * log(1 + (y-mu)^2 / (exp(f)*sigma^2)))
    }

    offset <- function(y,w)
        log(df)
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



StudentTLSS <- function(mu = 0, sigma = 1, df = 4){
    if (sigma <= 0 || df <= 0)
        stop(sQuote("sigma"), " and ", sQuote("df"),
             " must be greater than zero")
    RET <- list(mu = StudentTMu(mu = mu, sigma = sigma, df = df),
                sigma = StudentTSigma(mu = mu, sigma = sigma, df = df),
                df = StudentTDf(mu = mu, sigma = sigma, df = df))
    class(RET) <- "families"
    return(RET)
}
