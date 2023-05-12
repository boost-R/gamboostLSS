################################################################################
### Family wrapper for gamlss families




################################################################################
## constructor

## a wrapper to as.families (for compatibility reasons
gamlss.Families <- function(...)
    as.families(...)

as.families <- function(fname = "NO", stabilization = c("none", "MAD", "L2"),
                        mu = NULL, sigma = NULL, nu = NULL, tau = NULL,
                        mu.link = NULL, sigma.link = NULL, nu.link = NULL, 
                        tau.link = NULL) {
    
    ## require gamlss.dist
    if (!requireNamespace("gamlss.dist", quietly = TRUE))
        stop("Please install package 'gamlss.dist' for using gamlss families.")
    
    if (is.function(fname))
        fname <- fname()
    
    if (inherits(fname, "gamlss.family"))
        fname <- fname$family[1]
    
    if (mode(fname) != "character" && mode(fname) != "name")
        fname <- as.character(substitute(fname))
    
    gamlss.fam <- try(gamlss.dist::gamlss.family(fname), silent = TRUE)
    if (inherits(gamlss.fam, "try-error")){
      suffix <- substr(fname, nchar(fname)-1, nchar(fname))
      distName <- substr(fname, 1, nchar(fname)-2)
      isCens <- c(as.integer(gregexpr("r",suffix)),
                  as.integer(gregexpr("l",suffix)),
                  as.integer(gregexpr("i",suffix)))
      typeCens <- which(isCens > 0)
      typesCens <- c("right", "left", "interval")
      if(length(typeCens)==1 && as.integer(gregexpr("c",suffix))>0 && is.list(try(gamlss.dist::gamlss.family(distName), silent = TRUE))){
        print(paste0(fname, " specifies no valid gamlss family. Try to define a ",
                     "censored family with the following command first\nand then,",
                     "run your code again:"))
        print(paste0("gen.cens(family = '", distName, "', name = 'cens', type = '", typesCens[typeCens], "')"))
      }
      stop(sQuote("fname"), " specifies no valid gamlss family")
    }
    
    stabilization <- check_stabilization(stabilization)
    
    npar <- gamlss.fam$nopar
    switch(npar, {
        ## 1 parameter
        fun <- gamlss1parMu(mu = mu, fname = fname, mu.link = mu.link)
        warning("For boosting one-parametric families,",
                " please use the mboost package.")
        if (stabilization != "none")
            warning("Stabilization is ignored for one-parametric families.")
    }, {
        ## 2 parameters
        fun <- gamlss2parFam(mu = mu, sigma = sigma, mu.link = mu.link, 
                             sigma.link = sigma.link,
                             stabilization = stabilization, fname = fname)
    }, {
        ## 3 parameters
        fun <- gamlss3parFam(mu = mu, sigma = sigma, nu = nu, 
                             mu.link = mu.link, sigma.link = sigma.link, 
                             nu.link = nu.link,
                             stabilization = stabilization, fname = fname)
    }, {
        ## 4 parameters
        fun <- gamlss4parFam(mu = mu, sigma = sigma, nu = nu, tau = tau,
                             mu.link = mu.link, sigma.link = sigma.link,
                             nu.link = nu.link, tau.link = tau.link,
                             stabilization = stabilization, fname = fname)
    })
    fun
}


################################################################################
## 1 parameter

gamlss1parMu <- function(mu = NULL, fname = "EXP", mu.link = mu.link) {
    
    ## check if a specific link was provided 
    if (!is.null(mu.link)) fname_link <- eval(parse(text = paste(fname, 
                                                                "(mu.link = '", mu.link, "')", sep ="")))
    else fname_link <- fname
    
    FAM <- gamlss.dist::as.gamlss.family(fname_link)
    NAMEofFAMILY <- FAM$family
    dfun <- paste("gamlss.dist::d", fname, sep = "")
    pdf <- eval(parse(text = dfun))
    is.bdfamily  <- "bd" %in% names(formals(pdf))
    
    
    ## get the loss
    loss <- function(y, f, w = 1) {
        if (is.bdfamily) {
            bd <- rowSums(y)
            y <- y[,1]
            return( -pdf(x = y, mu = FAM$mu.linkinv(f), log = TRUE, bd = bd))
        }
        -pdf(x = y, mu = FAM$mu.linkinv(f), log = TRUE)
    }
    
    ## compute the risk
    risk <- function(y, f, w = 1) {
        sum(w * loss(y = y, f = f))
    }
    ## get the ngradient: mu is linkinv(f)
    ## we need dl/deta = dl/dmu*dmu/deta
    ngradient <- function(y, f, w = 1) {
        if (is.bdfamily) {
            if (!is.matrix(y)) stop("y should be a matrix for this family")
            bd <- rowSums(y)
            y <- y[,1]
            ngr <- FAM$dldm(y = y, mu = FAM$mu.linkinv(f), bd = bd) * FAM$mu.dr(eta = f)
        } else {
            ngr <- FAM$dldm(y = y, mu = FAM$mu.linkinv(f)) * FAM$mu.dr(eta = f)  
        }
        return(ngr)
    }
    
    ## get the offset -> we take the starting values of gamlss
    offset <- function(y, w = 1) {
        if (!is.null(mu)) {
            RET <- FAM$mu.linkfun(mu)
        } else {
            if (is.bdfamily) {
                if (!is.matrix(y)) stop("y should be a matrix for this family")
                bd <- rowSums(y)
                y <- y[,1]
            }
            eval(FAM$mu.initial)
            RET <- FAM$mu.linkfun(mean(mu))
        }
        return(RET)
    }
    
    mboost::Family(ngradient = ngradient, risk = risk, loss = loss,
                   response = function(f) FAM$mu.linkinv(f), offset = offset,
                   name = paste(FAM$family[2], "(mboost family)"))
}


################################################################################
## 2 parameters

gamlss2parMu <- function(mu = NULL, sigma = NULL, mu.link = NULL, FAM = FAM, 
                         stabilization = stabilization, fname = "NO") {
    
    NAMEofFAMILY <- FAM$family
    pdf <- get_pdf(fname)
    is.bdfamily  <- "bd" %in% names(formals(pdf))
    
    ## get the loss
    loss <- function(y, f, sigma,  w = 1) {
        if (is.bdfamily) {
            #if (!is.matrix(y)) stop("y should be a matrix for this family")
            bd <- rowSums(y)
            y <- y[,1]
            return( -pdf(x = y, mu = FAM$mu.linkinv(f), sigma = sigma, log = TRUE, bd = bd))
        }
        
        -pdf(x = y, mu = FAM$mu.linkinv(f), sigma = sigma, log = TRUE)
    }
    
    ## compute the risk
    risk <- function(y, f, w = 1) {
        sum(w * loss(y = y, f = f, sigma = sigma))
    }
    
    ## get the ngradient: mu is linkinv(f)
    ## we need dl/deta = dl/dmu*dmu/deta
    ngradient <- function(y, f, w = 1) {
        if (is.bdfamily) {
            if (!is.matrix(y)) stop("y should be a matrix for this family")
            bd <- rowSums(y)
            y <- y[,1]
            ngr <-  FAM$dldm(y = y, mu = FAM$mu.linkinv(f), sigma = sigma, bd = bd) * FAM$mu.dr(eta = f)
        } else {
            ngr <-  FAM$dldm(y = y, mu = FAM$mu.linkinv(f), sigma = sigma) * FAM$mu.dr(eta = f)  
        }
        ngr <- stabilize_ngradient(ngr, w = w, stabilization)
        ngr
    }
    
    ## get the offset -> we take the starting values of gamlss
    offset <- function(y, w = 1) {
        if (!is.null(mu)) {
            RET <- FAM$mu.linkfun(mu)
        } else {
            if (is.bdfamily) {
                if (!is.matrix(y)) stop("y should be a matrix for this family")
                bd <- rowSums(y)
                y <- y[,1]
            }
            eval(FAM$mu.initial)
            RET <- FAM$mu.linkfun(mean(mu))
        }
        return(RET)
    }
    
    Family(ngradient = ngradient, risk = risk, loss = loss,
           response = function(f) FAM$mu.linkinv(f), offset = offset,
           name = paste(FAM$family[2], "1st parameter (mu)"))
}


gamlss2parSigma <- function(mu = NULL, sigma = NULL,
                            stabilization = stabilization, fname = "NO", FAM = FAM) {
    
    NAMEofFAMILY <- FAM$family
    pdf <- get_pdf(fname)
    is.bdfamily  <- "bd" %in% names(formals(pdf))
    
    ## get the loss
    loss <- function(y, f, w = 1, mu ) {
        # check if bd exists in this family
        if (is.bdfamily) {
            bd <- rowSums(y)
            y <- y[,1]
            return( -pdf(x = y, mu = mu, sigma = FAM$sigma.linkinv(f), log = TRUE, bd = bd))
        }
        
        -pdf(x = y, mu = mu, sigma = FAM$sigma.linkinv(f), log = TRUE)
    }
    
    ## compute the risk
    risk <- function(y, f, w = 1) {
        sum(w * loss(y = y, f = f, mu = mu))
    }
    ## get the ngradient: sigma is linkinv(f)
    ## we need dl/deta = dl/dsigma*dsigma/deta
    ngradient <- function(y, f, w = 1) {
        
        # check if bd exists in this family
        if (is.bdfamily) {
            if (!is.matrix(y)) stop("y should be a matrix for this family")
            bd <- rowSums(y)
            y <- y[,1]
            ngr <- FAM$dldd(y = y, mu = mu, sigma = FAM$sigma.linkinv(f), bd = bd) * FAM$sigma.dr(eta = f)
            
        } else { 
	    ngr <- FAM$dldd(y = y, mu = mu, sigma = FAM$sigma.linkinv(f)) * FAM$sigma.dr(eta = f)
        }
        
        ngr <- stabilize_ngradient(ngr, w = w, stabilization)
        ngr
    }
    ## get the offset
    offset <- function(y, w = 1) {
        if (!is.null(sigma)) {
            RET <- FAM$sigma.linkfun(sigma)
        } else {
            eval(FAM$sigma.initial)
            RET <- FAM$sigma.linkfun(mean(sigma))
        }
        return(RET)
    }
    Family(ngradient = ngradient, risk = risk, loss = loss,
           response = function(f) FAM$sigma.linkinv(f), offset = offset,
           name = paste(FAM$family[2], "2nd parameter (sigma)"))
}

## Build the Families object
gamlss2parFam <- function(mu = NULL, sigma = NULL, mu.link = mu.link,  
                          sigma.link = sigma.link, stabilization, fname = "NO") {
    
    FAM <- gamlss.dist::as.gamlss.family(fname)
    # check if any link function was set    
    if (any(!is.null(mu.link), !is.null(sigma.link))) {
        # if some are null, set default
        if (is.null(mu.link)) mu.link <- FAM$mu.link 
        if (is.null(sigma.link)) sigma.link <- FAM$sigma.link 
        fname_link <- paste("gamlss.dist::",fname, "(mu.link = '", mu.link, "', ", 
                            "sigma.link = '", sigma.link, "'",  ")", sep ="")
        # build family with link
        FAM <- gamlss.dist::as.gamlss.family(eval(parse(text = fname_link)))
    } 
    Families(mu = gamlss2parMu(mu = mu, sigma = sigma, FAM = FAM, 
                               stabilization = stabilization, fname = fname),
             sigma = gamlss2parSigma(mu = mu, sigma = sigma, FAM = FAM, 
                                     stabilization = stabilization, fname = fname),
             qfun = get_qfun(fname),
             name = fname)
}

################################################################################
## 3 parameters

## sub-family for Mu
gamlss3parMu <- function(mu = NULL, sigma = NULL, nu = NULL, FAM = FAM,
                         stabilization = stabilization, fname = "TF") {
    
    
    NAMEofFAMILY <- FAM$family
    pdf <- get_pdf(fname)
    
    ## get the loss
    loss <- function(y, f, sigma, nu, w = 1) {
        -pdf(x = y, mu = FAM$mu.linkinv(f), sigma = sigma, nu = nu, log = TRUE)
    }
    
    ## compute the risk
    risk <- function(y, f, w = 1) {
        sum(w * loss(y = y, f = f, sigma = sigma, nu = nu))
    }
    ## get the ngradient: mu is linkinv(f)
    ## we need dl/deta = dl/dmu*dmu/deta
    ngradient <- function(y, f, w = 1) {
        if (FAM$type == "Mixed") {
            ngr <- FAM$dldm(y = y, mu = FAM$mu.linkinv(f), sigma = sigma) * FAM$mu.dr(eta = f)
        } else {
            ngr <- FAM$dldm(y = y, mu = FAM$mu.linkinv(f), sigma = sigma, nu = nu) * FAM$mu.dr(eta = f)
        }
        ngr <- stabilize_ngradient(ngr, w = w, stabilization)
        ngr
    }
    
    ## get the offset -> we take the starting values of gamlss
    offset <- function(y, w = 1) {
        if (!is.null(mu)) {
            RET <- FAM$mu.linkfun(mu)
        } else {
            eval(FAM$mu.initial)
            RET <- FAM$mu.linkfun(mean(mu))
        }
        return(RET)
    }
    Family(ngradient = ngradient, risk = risk, loss = loss,
           response = function(f) FAM$mu.linkinv(f), offset = offset,
           name = paste(FAM$family[2], "1st parameter (mu)"))
}


gamlss3parSigma <- function(mu = NULL, sigma = NULL, nu = NULL, FAM = FAM,
                            stabilization = stabilization, fname = "TF") {
    
    NAMEofFAMILY <- FAM$family
    pdf <- get_pdf(fname)
    
    ## get the loss
    loss <- function(y, f, w = 1, mu, nu) {
        -pdf(x = y, mu = mu, sigma = FAM$sigma.linkinv(f), nu = nu, log = TRUE)
    }
    
    ## compute the risk
    risk <- function(y, f, w = 1) {
        sum(w * loss(y = y, f = f, mu = mu, nu = nu))
    }
    ## get the ngradient: sigma is linkinv(f)
    ## we need dl/deta = dl/dsigma*dsigma/deta
    ngradient <- function(y, f, w = 1) {
        if (FAM$type == "Mixed") {
            ngr <- FAM$dldd(y = y, mu = mu, sigma = FAM$sigma.linkinv(f))  * FAM$sigma.dr(eta = f)
        } else {
            ngr <- FAM$dldd(y = y, mu = mu, sigma = FAM$sigma.linkinv(f), nu = nu) * FAM$sigma.dr(eta = f)
        }
        ngr <- stabilize_ngradient(ngr, w = w, stabilization)
        ngr
    }
    ## get the offset
    offset <- function(y, w = 1) {
        if (!is.null(sigma)) {
            RET <- FAM$sigma.linkfun(sigma)
        } else {
            eval(FAM$sigma.initial)
            RET <- FAM$sigma.linkfun(mean(sigma))
        }
        return(RET)
    }
    Family(ngradient = ngradient, risk = risk, loss = loss,
           response = function(f) FAM$sigma.linkinv(f), offset = offset,
           name = paste(FAM$family[2], "2nd parameter (sigma)"))
}

gamlss3parNu <- function(mu = NULL, sigma = NULL, nu = NULL,
                         stabilization = stabilization, fname = "TF", FAM = FAM) {
    
    NAMEofFAMILY <- FAM$family
    pdf <- get_pdf(fname)
    
    ## get the loss
    loss <- function(y, f, w = 1, mu, sigma) {
        -pdf(x = y, mu = mu, sigma = sigma, nu = FAM$nu.linkinv(f), log = TRUE)
    }
    
    ## compute the risk
    risk <- function(y, f, w = 1) {
        sum(w * loss(y = y, f = f, mu = mu, sigma = sigma))
    }
    ## get the ngradient: sigma is linkinv(f)
    ## we need dl/deta = dl/dsigma*dsigma/deta
    ngradient <- function(y, f, w = 1) {
        if (FAM$type == "Mixed") {
            ngr <- FAM$dldv(y = y, nu = FAM$nu.linkinv(f)) * FAM$nu.dr(eta = f)
        } else {
            ngr <- FAM$dldv(y = y, mu = mu, sigma = sigma, nu = FAM$nu.linkinv(f)) * FAM$nu.dr(eta = f)
        }
        ngr <- stabilize_ngradient(ngr, w = w, stabilization)
        ngr
    }
    ## get the offset
    offset <- function(y, w = 1) {
        if (!is.null(nu)) {
            RET <- FAM$nu.linkfun(nu)
        } else {
            eval(FAM$nu.initial)
            RET <- FAM$nu.linkfun(mean(nu))
        }
        return(RET)
    }
    Family(ngradient = ngradient, risk = risk, loss = loss,
           response = function(f) FAM$nu.linkinv(f),
           offset = offset, name = paste(FAM$family[2], "3rd parameter (nu)"))
}

## Build the Families object
gamlss3parFam <- function(mu = NULL, sigma = NULL, nu = NULL, 
                          mu.link = NULL, sigma.link = NULL, nu.link = NULL, 
                          stabilization = stabilization, fname = "TF") {
    
    FAM <- gamlss.dist::as.gamlss.family(fname)
    # check if any link function was set    
    if(any(!is.null(mu.link), !is.null(sigma.link), !is.null(nu.link))){
        # if some are null, set default
        if (is.null(mu.link)) mu.link <- FAM$mu.link 
        if (is.null(sigma.link)) sigma.link <- FAM$sigma.link 
        if (is.null(nu.link)) nu.link <- FAM$nu.link 
        
        fname_link <- paste("gamlss.dist::",fname, "(mu.link = '", mu.link, "', ", 
                            "sigma.link = '", sigma.link, "',",  
                            "nu.link = '", nu.link, "')", sep ="")
        # build family with link
        FAM <- gamlss.dist::as.gamlss.family(eval(parse(text = fname_link)))
    }   
    
    Families(mu = gamlss3parMu(mu = mu, sigma = sigma, nu = nu, 
                               stabilization = stabilization, fname = fname, FAM = FAM),
             sigma = gamlss3parSigma(mu = mu, sigma = sigma, nu = nu, 
                                     stabilization = stabilization, fname = fname, FAM = FAM),
             nu = gamlss3parNu(mu = mu, sigma = sigma, nu = nu, 
                               stabilization = stabilization, fname = fname, FAM = FAM),
             qfun = get_qfun(fname),
             name = fname)
}


################################################################################
## 4 parameters

gamlss4parMu <- function(mu = NULL, sigma = NULL, nu = NULL, tau = NULL,
                         stabilization = stabilization, fname = "BCPE", FAM = FAM) {
    
    NAMEofFAMILY <- FAM$family
    pdf <- get_pdf(fname)
    
    ## get the loss
    loss <- function(y, f, sigma, nu, tau, w = 1) {
        -pdf(x = y, mu = FAM$mu.linkinv(f), sigma = sigma, nu = nu, tau = tau, log = TRUE)
    }
    ## compute the risk
    risk <- function(y, f, w = 1) {
        ## get the ngradient: mu is linkinv(f)
        sum(w * loss(y = y, f = f, sigma = sigma, nu = nu, tau = tau))
    }
    ## we need dl/deta = dl/dmu*dmu/deta
    ngradient <- function(y, f, w = 1) {
        if (FAM$type == "Mixed") {
            ngr <- FAM$dldm(y = y, mu = FAM$mu.linkinv(f), sigma = sigma) * FAM$mu.dr(eta = f)
        } else {
            ngr <- FAM$dldm(y = y, mu = FAM$mu.linkinv(f), sigma = sigma, nu = nu, tau = tau) *
                FAM$mu.dr(eta = f)
        }
        ngr <- stabilize_ngradient(ngr, w = w, stabilization)
        ngr
    }
    ## get the offset -> we take the starting values of gamlss
    offset <- function(y, w = 1) {
        if (!is.null(mu)) {
            RET <- FAM$mu.linkfun(mu)
        } else {
            eval(FAM$mu.initial)
            RET <- FAM$mu.linkfun(mean(mu))
        }
        return(RET)
    }
    Family(ngradient = ngradient, risk = risk, loss = loss,
           response = function(f) FAM$mu.linkinv(f), offset = offset,
           name = paste(FAM$family[2], "1st parameter (mu)"))
}


gamlss4parSigma <- function(mu = NULL, sigma = NULL, nu = NULL, tau = NULL, FAM = FAM,
                            stabilization = stabilization, fname = "BCPE") {
    
    NAMEofFAMILY <- FAM$family
    pdf <- get_pdf(fname)
    
    ## get the loss
    loss <- function(y, f, w = 1, mu, nu, tau) {
        -pdf(x = y, mu = mu, sigma = FAM$sigma.linkinv(f), nu = nu, tau = tau,
             log = TRUE)
    }
    
    ## compute the risk
    risk <- function(y, f, w = 1) {
        sum(w * loss(y = y, f = f, mu = mu, nu = nu, tau = tau))
    }
    ## get the ngradient: sigma is linkinv(f)
    ## we need dl/deta = dl/dsigma*dsigma/deta
    ngradient <- function(y, f, w = 1) {
        if (FAM$type == "Mixed") {
            ngr <- FAM$dldm(y = y, mu = mu, sigma = FAM$sigma.linkinv(f)) * FAM$sigma.dr(eta = f)
        } else {
            ngr <-  FAM$dldd(y = y, mu = mu, sigma = FAM$sigma.linkinv(f), nu = nu,
                             tau = tau) * FAM$sigma.dr(eta = f)
        }
        ngr <- stabilize_ngradient(ngr, w = w, stabilization)
        ngr
    }
    ## get the offset
    offset <- function(y, w = 1) {
        if (!is.null(sigma)) {
            RET <- FAM$sigma.linkfun(sigma)
        } else {
            eval(FAM$sigma.initial)
            RET <- FAM$sigma.linkfun(mean(sigma))
        }
        return(RET)
    }
    Family(ngradient = ngradient, risk = risk, loss = loss,
           response = function(f) FAM$sigma.linkinv(f), offset = offset,
           name = paste(FAM$family[2], "2nd parameter (sigma)"))
}

gamlss4parNu <- function(mu = NULL, sigma = NULL, nu = NULL, tau = NULL, FAM = FAM,
                         stabilization = stabilization, fname = "BCPE") {
    
    NAMEofFAMILY <- FAM$family
    pdf <- get_pdf(fname)
    
    ## get the loss
    loss <- function(y, f, w = 1, mu, sigma, tau) {
        -pdf(x = y, mu = mu, sigma = sigma, nu = FAM$nu.linkinv(f), tau = tau,
             log = TRUE)
    }
    
    ## compute the risk
    risk <- function(y, f, w = 1) {
        sum(w * loss(y = y, f = f, mu = mu, sigma = sigma, tau = tau))
    }
    ## get the ngradient: sigma is linkinv(f)
    ## we need dl/deta = dl/dsigma*dsigma/deta
    ngradient <- function(y, f, w = 1) {
        if (FAM$type == "Mixed") {
            ngr <- FAM$dldv(y = y, nu = FAM$nu.linkinv(f), tau = tau) * FAM$nu.dr(eta = f)
        } else {
            ngr <- FAM$dldv(y = y, mu = mu, sigma = sigma, nu = FAM$nu.linkinv(f),
                            tau = tau) * FAM$nu.dr(eta = f)
        }
        ngr <- stabilize_ngradient(ngr, w = w, stabilization)
        ngr
    }
    ## get the offset
    offset <- function(y, w = 1) {
        if (!is.null(nu)) {
            RET <- FAM$nu.linkfun(nu)
        } else {
            eval(FAM$nu.initial)
            RET <- FAM$nu.linkfun(mean(nu))
        }
        return(RET)
    }
    Family(ngradient = ngradient, risk = risk, loss = loss,
           response = function(f) FAM$nu.linkinv(f), offset = offset,
           name = paste(FAM$family[2], "3rd parameter (nu)"))
}



gamlss4parTau <- function(mu = NULL, sigma = NULL, nu = NULL, tau = NULL, FAM = FAM,
                          stabilization = stabilization, fname = "BCPE") {
    
    NAMEofFAMILY <- FAM$family
    pdf <- get_pdf(fname)
    
    ## get the loss
    loss <- function(y, f, w = 1, mu, sigma, nu) {
        -pdf(x = y, mu = mu, sigma = sigma, nu = nu, tau = FAM$tau.linkinv(f),
             log = TRUE)
    }
    
    ## compute the risk
    risk <- function(y, f, w = 1) {
        sum(w * loss(y = y, f = f, mu = mu, sigma = sigma, nu = nu))
    }
    ## get the ngradient
    ngradient <- function(y, f, w = 1) {
        if (FAM$type == "Mixed") {
            ngr <- FAM$dldt(y = y, nu = nu, tau =  FAM$tau.linkinv(f)) * FAM$tau.dr(eta = f)
        } else {
            ngr <- FAM$dldt(y = y, mu = mu, sigma = sigma, tau = FAM$tau.linkinv(f),
                            nu = nu) * FAM$tau.dr(eta = f)
        }
        ngr <- stabilize_ngradient(ngr, w = w, stabilization)
        ngr
    }
    ## get the offset
    offset <- function(y, w = 1) {
        if (!is.null(tau)) {
            RET <- FAM$tau.linkfun(tau)
        } else {
            eval(FAM$tau.initial)
            RET <- FAM$tau.linkfun(mean(tau))
        }
        return(RET)
    }
    Family(ngradient = ngradient, risk = risk, loss = loss,
           response = function(f) FAM$tau.linkinv(f), offset = offset,
           name = paste(FAM$family[2], "4th parameter (tau)"))
}

## Build the Families object
gamlss4parFam <- function(mu = NULL, sigma = NULL, nu = NULL, tau = NULL, 
                          mu.link = NULL, sigma.link = NULL, nu.link = NULL, tau.link = NULL, 
                          stabilization = stabilization, fname = "BCPE") {
    
    FAM <- gamlss.dist::as.gamlss.family(fname)
    # check if any link function was set    
    if (any(!is.null(mu.link), !is.null(sigma.link), !is.null(nu.link), 
           !is.null(tau.link))) {
        # if some are null, set default
        if (is.null(mu.link)) mu.link <- FAM$mu.link 
        if (is.null(sigma.link)) sigma.link <- FAM$sigma.link 
        if (is.null(nu.link)) nu.link <- FAM$nu.link 
        if (is.null(tau.link)) tau.link <- FAM$tau.link 
        
        
        fname_link <- paste("gamlss.dist::",fname, "(mu.link = '", mu.link, "', ", 
                            "sigma.link = '", sigma.link, "',",  
                            "nu.link = '", nu.link, "',",
                            "tau.link = '", tau.link, "')", sep ="")
        # build family with link
        FAM <- gamlss.dist::as.gamlss.family(eval(parse(text = fname_link)))
    }   
    
    
    Families(mu = gamlss4parMu(mu = mu, sigma = sigma, nu = nu, tau = tau, 
                               stabilization = stabilization, fname = fname, FAM = FAM),
             sigma = gamlss4parSigma(mu = mu, sigma = sigma, nu = nu, tau = tau, 
                                     stabilization = stabilization, fname = fname, FAM = FAM),
             nu = gamlss4parNu(mu = mu, sigma = sigma, nu = nu, tau = tau, 
                               stabilization = stabilization, fname = fname, FAM = FAM),
             tau = gamlss4parTau(mu = mu, sigma = sigma, nu = nu, tau = tau, 
                                 stabilization = stabilization, fname = fname, FAM = FAM),
             qfun = get_qfun(fname),
             name = fname)
}
