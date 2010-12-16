###
# Constructor Function

Families <- function(...){
    RET <- list(...)
    class(RET) <- "families"
    ## check if response function is not specified
    check <- sapply(RET, function(x){
        bdy <- body(x@response)
        (length(x@response(c(-1,0,1))) != 3 && class(bdy) != "call" &&
         length(bdy) == 1 && is.na(bdy))
    })
    if (any(check))
        stop("response function not specified in families for:\n\t",
             paste(names(RET)[check], collapse =", "))

    return(RET)
}

###
# Negative Binomial LSS Family

NBinomialMu <- function(mu = NULL, sigma = NULL) {
    loss <- function(sigma, y, f)
        -(lgamma(y + sigma) - lgamma(sigma) - lgamma(y + 1) + sigma * log(sigma)
          + y * f - (y + sigma) * log(exp(f) + sigma))
    risk <- function(y, f, w = 1){
        RET <- sum(w * loss(y = y, f = f, sigma = sigma))
        return(RET)
    }
    ngradient <- function(y, f, w = 1) {
        RET <- y - (y + sigma)/(exp(f) + sigma) * exp(f)
        return(RET)
    }
    offset <- function(y, w){
        if (!is.null(mu)){
            RET <- log(mu)
        } else {
            if (is.null(sigma))
                sigma <<- mean(y)^2 / (sd(y) - mean(y))
            ### look for starting value of f = log(mu) in "interval"
            ### i.e. mu possibly ranges from 1e-10 to 1e10
            RET <- optimize(risk, interval = c(log(1e-10), log(1e10)), y = y, w = w)$minimum
        }
        return(RET)
    }

    Family(ngradient = ngradient, offset = offset, risk = risk, loss = loss,
           response = function(f) exp(f),
           check_y = function(y) {
               stopifnot(all.equal(unique(y - floor(y)), 0))
               ## <FIXME> Does not check if y >= 0
               y
           }, name = "Negative Negative-Binomial Likelihood: mu (log link)")
}

NBinomialSigma <- function(mu = NULL, sigma = NULL) {
    # this family boosts sigma therefore f is sigma
    loss <- function(mu, y, f)
        -(lgamma(y + exp(f)) - lgamma(exp(f)) - lgamma(y + 1) + exp(f) * f +
          y * log(mu) - (y + exp(f)) * log(mu + exp(f)))
    risk <- function(y, f, w = 1){
        RET <- sum(w * loss(y = y, f = f, mu = mu))
        return(RET)
    }
    ngradient <- function(y, f, w = 1) {       # f is sigma !
        RET <- exp(f)*(digamma(y +exp(f)) - digamma(exp(f)) + log(exp(f)) + 1 -
                       log(mu +exp(f)) - (exp(f) + y)/(mu +exp(f)))
        return(RET)
    }
    offset <- function(y,w){
        if (!is.null(sigma)){
            RET <- log(sigma)
        } else {
            if (is.null(mu))
                mu <<- mean(y)
            ### look for starting value of f = log(sigma) in "interval"
            ### i.e. sigma possibly ranges from 1e-10 to 1e10
            RET <- optimize(risk, interval = c(log(1e-10), log(1e10)), y = y, w = w)$minimum
        }
        return(RET)
    }

    Family(ngradient = ngradient, offset = offset, risk = risk, loss = loss,
           response = function(f) exp(f),
           check_y = function(y) {
               stopifnot(all.equal(unique(y - floor(y)), 0))
               ## <FIXME> Does not check if y >= 0
               y
           }, name = "Negative Negative-Binomial Likelihood: sigma (log link)")
}


NBinomialLSS <- function(mu = NULL, sigma = NULL){
    if ((!is.null(sigma) && sigma <= 0) || (!is.null(mu) && mu <= 0))
        stop(sQuote("sigma"), " and ", sQuote("mu"),
             " must be greater than zero")
    Families(mu = NBinomialMu(mu = mu, sigma = sigma),
             sigma = NBinomialSigma(mu = mu, sigma = sigma))
}



###
# T-distribution LSS Family

StudentTMu <- function(mu = NULL, sigma = NULL, df = NULL) {
    loss <- function(df, sigma,y, f){
        -1 * (lgamma((df+1)/2) - log(sigma) - lgamma(1/2) - lgamma(df/2) - 0.5 *
              log(df) - (df+1)/2 * log(1 + (y-f)^2 / (df * sigma^2)))
    }
    risk <- function(y, f, w = 1){
        sum(w * loss(y = y, f = f, df = df, sigma = sigma))
    }
    ngradient <- function(y, f, w = 1) {
        (df+1)*(y-f)/(df*sigma^2 +(y-f)^2)
     }

    offset <- function(y, w){
        if (!is.null(mu)){
            RET <- mu
        } else {
            if (is.null(sigma))
                sigma <<- 1
            if (is.null(df))
                df <<- 4
            RET <- optimize(risk, interval = range(y), y = y, w = w)$minimum
        }
        return(RET)
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
    risk <- function(y, f, w = 1){
        sum(w * loss(y = y, f = f, df = df, mu = mu))
    }
    ngradient <- function(y, f, w = 1) {
        (-1 + (df+1)/(df*exp(2*f)/(y-mu)^2 + 1))
    }
    offset <- function(y, w){
        if (!is.null(sigma)){
            RET <- log(sigma)
        } else {
            if (is.null(mu))
                mu <<- mean(y)
            if (is.null(df))
                df <<- 4
            ### look for starting value of f = log(sigma) in "interval"
            ### i.e. sigma possibly ranges from 1e-10 to 1e10
            RET <- optimize(risk, interval = c(log(1e-10), log(1e10)), y = y, w = w)$minimum
        }
        return(RET)
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
    risk <- function(y, f, w = 1){
        sum(w * loss(y = y, f = f, sigma = sigma, mu = mu))
    }
    ngradient <- function(y, f, w = 1) {
        exp(f)/2 * (digamma((exp(f)+1)/2)-digamma(exp(f)/2)) - 0.5 -
            (exp(f)/2 * log(1+ (y-mu)^2/(exp(f)*sigma^2)) -
             (y-mu)^2/(sigma^2 + (y-mu)^2/exp(f)) * (exp(-f) +1)/2 )
    }
    offset <- function(y, w){
        if (!is.null(df)){
            RET <- log(df)
        } else {
            if (is.null(mu))
                mu <<- mean(y)
            if (is.null(sigma))
                sigma <<- 1
            ### look for starting value of f = log(df) in "interval"
            ### i.e. df possibly ranges from 1e-10 to 1e10
            RET <- optimize(risk, interval = c(log(1e-10), log(1e10)), y = y, w = w)$minimum
        }
        return(RET)
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



StudentTLSS <- function(mu = NULL, sigma = NULL, df = NULL){
    if ((!is.null(sigma) && sigma <= 0) || (!is.null(df) && df <= 0))
        stop(sQuote("sigma"), " and ", sQuote("df"),
             " must be greater than zero")
    Families(mu = StudentTMu(mu = mu, sigma = sigma, df = df),
             sigma = StudentTSigma(mu = mu, sigma = sigma, df = df),
             df = StudentTDf(mu = mu, sigma = sigma, df = df))
}



###
# Log-Normal LSS Family

LogNormalMu <- function (mu = NULL, sigma = NULL){
    loss <- function(sigma, y, f) {
        logfw <- function(pred)
            dnorm(pred, log = TRUE)
        logSw <- function(pred)
            pnorm(pred, lower.tail = FALSE, log.p = TRUE)
        eta <- (log(y[,1]) - f)/sigma
        -y[,2] * (logfw(eta) - log(sigma)) - (1 - y[,2]) * logSw(eta)
    }
    risk <- function(y, f, w = 1)
        sum(w * loss(y = y, f = f, sigma = sigma))
    ngradient <- function(y, f, w = 1) {
        eta <- (log(y[,1]) - f)/sigma
        (y[,2] * eta + (1 - y[,2]) * dnorm(eta)/(1 - pnorm(eta)))/sigma
    }
    offset <- function(y, w){
        if (!is.null(mu)){
            RET <- mu
        } else {
            if (is.null(sigma))
                sigma <<- 1
            RET <- optimize(risk, interval = c(log(1e-10), log(1e10)),
                            y = y, w = w)$minimum
        }
        return(RET)
    }

    Family(ngradient = ngradient, risk = risk, offset = offset, loss = loss,
           response = function(f) f,
           check_y = function(y) {
               if (!inherits(y, "Surv"))
                   stop("response is not an object of class ", sQuote("Surv"),
                        " but ", sQuote("family = Lognormal()"))
               y
           }, name = "Log-Normal AFT Model: mu (id link)")
}

LogNormalSigma <- function(mu = NULL, sigma = NULL){
    loss <- function(mu, y, f) {
        logfw <- function(pred)
            dnorm(pred, log = TRUE)
        logSw <- function(pred)
            pnorm(pred, lower.tail = FALSE, log.p = TRUE)
        eta <- (log(y[,1]) - mu) / exp(f)
        -y[,2] * (logfw(eta) - f) - (1 - y[,2]) * logSw(eta)
    }
    risk <- function(y, f, w = 1)
        sum(w * loss(y = y, f = f, mu = mu))
    ngradient <- function(y, f, w = 1) {
        eta <- (log(y[,1]) - mu)/exp(f)
        -(y[,2] - y[,2]*eta^2 + (y[,2]-1)*eta*dnorm(eta)/(1-pnorm(eta)))
    }
    offset <- function(y, w){
        if (!is.null(sigma)){
            RET <- log(sigma)
        } else {
            if (is.null(mu))
                mu <<- 0
            ## look for starting value of f = log(sigma) in "interval"
            ## i.e. sigma possibly ranges from 1e-10 to 1e10
            RET <- optimize(risk, interval = c(log(1e-10), log(1e10)),
                            y = y, w = w)$minimum
        }
        return(RET)
    }

    Family(ngradient = ngradient, risk = risk, offset = offset, loss = loss,
           response = function(f) exp(f),
           check_y = function(y) {
               if (!inherits(y, "Surv"))
                   stop("response is not an object of class ", sQuote("Surv"),
                        " but ", sQuote("family = Lognormal()"))
               y
           }, name = "Log-Normal AFT Model: sigma (log link)")
}

LogNormalLSS <- function(mu = NULL, sigma = NULL){
    if ((!is.null(sigma) && sigma <= 0))
        stop(sQuote("sigma"), " must be greater than zero")
    Families(mu = LogNormalMu(mu = mu, sigma = sigma),
             sigma = LogNormalSigma(mu = mu, sigma = sigma))
}


###
# LogLog LSS Family

LogLogMu <- function (mu = NULL, sigma = NULL){
    loss <- function(sigma, y, f) {
        logfw <- function(pred)
            dlogis(pred, log = TRUE)
            #pred - 2 * (1 + exp(pred))
        logSw <- function(pred)
            plogis(pred, lower.tail = FALSE, log.p = TRUE)
            #1/(1 + exp(pred))
        eta <- (log(y[,1]) - f)/sigma
        -y[,2] * (logfw(eta) - log(sigma)) - (1 - y[,2]) * logSw(eta)
    }

    risk <- function(y, f, w = 1)
        sum(w * loss(y = y, f = f, sigma = sigma))

    ngradient <- function(y, f, w = 1) {
        eta <- (log(y[,1]) - f)/sigma
        nom <- (exp(-eta) + 1)
        (y[,2] * (2/nom - 1) + (1 - y[,2])/nom)/sigma
    }
    offset <- function(y, w){
        if (!is.null(mu)){
            RET <- mu
        } else {
            if (is.null(sigma))
                sigma <<- 1
            RET <- optimize(risk, interval = c(log(1e-10), log(1e10)),
                            y = y, w = w)$minimum
        }
        return(RET)
    }


    Family(ngradient = ngradient, risk = risk, offset = offset, loss = loss,
           response = function(f) f,
           check_y = function(y) {
               if (!inherits(y, "Surv"))
                   stop("response is not an object of class ", sQuote("Surv"),
                        " but ", sQuote("family = Loglog()"))
               y
           }, name = "Log-Logistic AFT Model: mu (id link)")
}

LogLogSigma <- function (mu = NULL, sigma = NULL){
    loss <- function(mu, y, f) {
        logfw <- function(pred)
            dlogis(pred, log = TRUE)
            #exp(pred)/(1 + exp(pred))^2
        logSw <- function(pred)
            pnorm(pred, lower.tail = FALSE, log.p = TRUE)
            #1/(1 + exp(pred))
        eta <- (log(y[,1]) - mu)/exp(f)
        -y[,2] * (logfw(eta) - f) - (1 - y[,2]) * logSw(eta)
    }
    risk <- function(y, f, w = 1)
        sum(w * loss(y = y, f = f, mu = mu))
    ngradient <- function(y, f, w = 1) {
        eta <- (log(y[,1]) - mu)/exp(f)
        -(y[,2] + y[,2]*eta -(y[,2]+1)*eta/(1+exp(-eta)))
    }
    offset <- function(y, w){
        if (!is.null(sigma)){
            RET <- log(sigma)
        } else {
            if (is.null(mu))
                mu <<- 0
            ## look for starting value of f = log(sigma) in "interval"
            ## i.e. sigma possibly ranges from 1e-10 to 1e10
            RET <- optimize(risk, interval = c(log(1e-10), log(1e10)), y = y, w = w)$minimum
        }
        return(RET)
    }

    Family(ngradient = ngradient, risk = risk, offset = offset, loss = loss,
           response = function(f) exp(f),
           check_y = function(y) {
               if (!inherits(y, "Surv"))
                   stop("response is not an object of class ", sQuote("Surv"),
                        " but ", sQuote("family = Loglog()"))
               y
           },  name = "Log-Logistic AFT Model: sigma (log link)")
}

LogLogLSS <- function(mu = NULL, sigma = NULL){
    if ((!is.null(sigma) && sigma <= 0))
        stop(sQuote("sigma"), " must be greater than zero")
    Families(mu = LogLogMu(mu = mu, sigma = sigma),
             sigma = LogLogSigma(mu = mu, sigma = sigma))
}


###
# Weibull LSS Family
WeibullMu <- function (mu = NULL, sigma = NULL){
    loss <- function(sigma, y, f) {
        logfw <- function(pred)
            pred - exp(pred)
        logSw <- function(pred)
            -exp(pred)
        eta <- (log(y[,1]) - f)/sigma
        -y[,2] * (logfw(eta) -log(sigma)) - (1 - y[,2]) * logSw(eta)
    }
    risk <- function(y, f, w = 1)
        sum(w * loss(y = y, f = f, sigma = sigma))
    ngradient <- function(y, f, w = 1){
        eta <- (log(y[,1]) - f)/sigma
        (y[,2] * (exp(eta) - 1) + (1 - y[,2]) * exp(eta))/sigma
    }
    offset <- function(y, w){
        if (!is.null(mu)){
            RET <- log(mu)
        } else {
            if (is.null(sigma))
                sigma <<- 1
            RET <- optimize(risk, interval = c(0, max(y[,1], na.rm = TRUE)),
                            y = y, w = w)$minimum
        }
        return(RET)
    }

    Family(ngradient = ngradient, risk = risk, offset = offset, loss = loss,
           response = function(f) f,
           check_y = function(y) {
               if (!inherits(y, "Surv"))
                   stop("response is not an object of class ", sQuote("Surv"),
                        " but ", sQuote("family = Weibull()"))
               y
           }, name = "Weibull AFT Model: mu (id link)")
}


WeibullSigma <- function (mu = NULL, sigma = NULL){
    loss <- function(mu, y, f) {
        logfw <- function(pred)
            pred - exp(pred)
        logSw <- function(pred)
            -exp(pred)
        eta <- (log(y[,1]) - mu)/exp(f)
        -y[,2] * (logfw(eta) - f) - (1 - y[,2]) * logSw(eta)
        }
    risk <- function(y, f, w = 1)
        sum(w * loss(y = y, f = f, mu = mu))
    ngradient <- function(y, f, w = 1) {
        eta <- (log(y[,1]) - mu)/exp(f)
        -(y[,2] * (eta + 1) - eta * exp(eta))
    }
    offset <- function(y, w){
        if (!is.null(sigma)){
            RET <- log(sigma)
        } else {
            if (is.null(mu))
                mu <<- 0
            ## look for starting value of f = log(sigma) in "interval"
            ## i.e. sigma possibly ranges from 1e-10 to 1e10
            RET <- optimize(risk, interval = c(log(1e-10), log(1e10)),
                            y = y, w = w)$minimum
        }
        return(RET)
    }

    Family(ngradient = ngradient, risk = risk, offset = offset, loss = loss,
           response = function(f) exp(f),
           check_y = function(y) {
               if (!inherits(y, "Surv"))
                   stop("response is not an object of class ", sQuote("Surv"),
                        " but ", sQuote("family = Weibull()"))
               y
           }, name = "Weibull AFT Model: sigma (log link)")
}

WeibullLSS <- function(mu = NULL, sigma = NULL){
    if ((!is.null(sigma) && sigma <= 0))
        stop(sQuote("sigma"), " must be greater than zero")
    Families(mu = WeibullMu(mu = mu, sigma = sigma),
             sigma = WeibullSigma(mu = mu, sigma = sigma))
}
