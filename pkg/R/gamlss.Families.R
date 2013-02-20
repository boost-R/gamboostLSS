
#
#
#
# Family wrapper for gamlss families



#-------------- 1 parameter 

gamlss1parMu <- function(mu = NULL, fname = "EXP"){
  
  FAM <- as.gamlss.family(fname)
  NAMEofFAMILY <- FAM$family
  dfun <- paste("d", fname, sep = "")
  pdf <- eval(parse(text = dfun))
  
  # get the loss
  loss <- function(y, f,  w = 1){
    - pdf(x = y, mu = FAM$mu.linkinv(f),  log = TRUE) 
  }
  
  # compute the risk
  risk <- function(y, f, w = 1) {
    sum(w * loss(y = y, f = f))
  }
  
  # get the ngradient:  mu is linkinv(f)
  # we need dl/deta = dl/dmu*dmu/deta 
  ngradient <- function(y, f, w = 1){
    FAM$dldm(y = y, mu = FAM$mu.linkinv(f))*FAM$mu.dr(eta = f)
  }
  
  # get the offset -> we take the starting values of gamlss
  offset <- function(y, w = 1)
  {
    if (!is.null(mu)){
      RET <- FAM$mu.linkfun(mu)
    } else {
      eval(FAM$mu.initial)
      RET <- FAM$mu.linkfun(mean(mu))
    }
    return(RET)
  }
  
  mboost:::Family(ngradient = ngradient, risk = risk, loss = loss,
         response = function(f) FAM$mu.linkinv(f), offset = offset, 
         name = paste(FAM$family[2], "mboost family"))
}


#-------------- 2 parameters 

gamlss2parMu <- function(mu = NULL,  sigma = NULL, fname="NO"){
  
  FAM <- as.gamlss.family(fname)
  NAMEofFAMILY <- FAM$family
  dfun <- paste("d", fname, sep = "")
  pdf <- eval(parse(text = dfun))
  
  # get the loss
  loss <- function(y, f, sigma,  w = 1){
    - pdf(x = y, mu = FAM$mu.linkinv(f), sigma = sigma, log = TRUE) 
  }
  
  # compute the risk
  risk <- function(y, f, w = 1) {
    sum(w * loss(y = y, f = f, sigma = sigma))
  }
  
  # get the ngradient:  mu is linkinv(f)
  # we need dl/deta = dl/dmu*dmu/deta 
  ngradient <- function(y, f, w = 1){
    FAM$dldm(y = y, mu = FAM$mu.linkinv(f), sigma = sigma)*FAM$mu.dr(eta = f)
  }  
  
  # get the offset -> we take the starting values of gamlss
  offset <- function(y, w = 1)
  {
    if (!is.null(mu)){
      RET <- FAM$mu.linkfun(mu)
    } else {
      eval(FAM$mu.initial)
      RET <- FAM$mu.linkfun(mean(mu))
    }
    return(RET)
  }
  mboost:::Family(ngradient = ngradient, risk = risk, loss = loss,
         response = function(f) FAM$mu.linkinv(f), offset=offset, 
         name = paste(FAM$family[2], "1st parameter (mu)"))
}


gamlss2parSigma <- function(mu = NULL,  sigma = NULL, fname="NO")
{
  FAM <- as.gamlss.family(fname)
  NAMEofFAMILY <- FAM$family
  dfun <- paste("d", fname, sep = "")
  pdf <- eval(parse(text = dfun))
  
  # get the loss
  loss <- function(y, f, w = 1, mu){
    - pdf(x = y, mu = mu, sigma = FAM$sigma.linkinv(f), log=TRUE) 
  }
  
  # compute the risk
  risk <- function(y, f, w = 1) {
    sum(w * loss(y = y, f = f, mu = mu))
  }
  # get the ngradient:  sigma is linkinv(f)
  # we need dl/deta = dl/dsigma*dsigma/deta
  ngradient <- function(y, f, w = 1){
    FAM$dldd(y = y, mu = mu, sigma = FAM$sigma.linkinv(f))*FAM$sigma.dr(eta = f)
  }
  # get the offset
  offset <- function(y, w = 1)
  {
    if (!is.null(sigma)){
      RET <- FAM$sigma.linkfun(sigma)
    } else {
      eval(FAM$sigma.initial)
      RET <- FAM$sigma.linkfun(mean(sigma))
    }
    return(RET)
  }
  mboost:::Family(ngradient = ngradient, risk = risk, loss = loss,
         response = function(f) FAM$sigma.linkinv(f), offset = offset, 
         name = paste(FAM$family[2], "2nd parameter (sigma)"))
}

# Build the Families object
gamlss2parFam <- function(mu = NULL, sigma = NULL, fname = "NO")
{ 
  Families(mu = gamlss2parMu(mu = mu, sigma = sigma, fname = fname), 
           sigma = gamlss2parSigma(mu = mu, sigma = sigma, fname = fname))
}

#----------------- 3 parameters 

# sub-family for Mu
gamlss3parMu <- function(mu = NULL,  sigma = NULL, nu = NULL, fname= "TF"){
  
  FAM <- as.gamlss.family(fname)
  NAMEofFAMILY <- FAM$family
  dfun <- paste("d", fname, sep = "")
  pdf <- eval(parse(text = dfun))
  
  # get the loss
  loss <- function(y, f, sigma, nu,  w = 1){
    - pdf(x = y, mu = FAM$mu.linkinv(f), sigma = sigma, nu = nu,  log=TRUE) 
  }
  
  # compute the risk
  risk <- function(y, f, w = 1) {
    sum(w * loss(y = y, f = f, sigma = sigma, nu = nu))
  }
  # get the ngradient:  mu is linkinv(f)
  # we need dl/deta = dl/dmu*dmu/deta 
  ngradient <- function(y, f, w = 1){
    FAM$dldm(y = y, mu = FAM$mu.linkinv(f), sigma = sigma, nu = nu)*FAM$mu.dr(eta = f)
  }
  # get the offset -> we take the starting values of gamlss
  offset <- function(y, w = 1)
  {
    if (!is.null(mu)){
      RET <- FAM$mu.linkfun(mu)
    } else {
      eval(FAM$mu.initial)
      RET <- FAM$mu.linkfun(mean(mu))
    }
    return(RET)
  }
  mboost:::Family(ngradient = ngradient, risk = risk, loss = loss,
         response = function(f) FAM$mu.linkinv(f), offset = offset, 
         name = paste(FAM$family[2], "1st parameter (mu)"))
}


gamlss3parSigma <- function(mu = NULL,  sigma = NULL, nu = NULL, fname = "TF")
{
  FAM <- as.gamlss.family(fname)
  NAMEofFAMILY <- FAM$family
  dfun <- paste("d", fname, sep = "")
  pdf <- eval(parse(text = dfun))
  
  # get the loss
  loss <- function(y, f, w = 1, mu, nu){
    - pdf(x = y, mu = mu, sigma = FAM$sigma.linkinv(f), nu = nu,  log = TRUE) 
  }
  
  # compute the risk
  risk <- function(y, f, w = 1) {
    sum(w * loss(y = y, f = f, mu = mu, nu = nu))
  }
  # get the ngradient:  sigma is linkinv(f)
  # we need dl/deta = dl/dsigma*dsigma/deta
  ngradient <- function(y, f, w = 1){
    FAM$dldd(y = y, mu = mu, sigma = FAM$sigma.linkinv(f), nu = nu)*FAM$sigma.dr(eta = f)
  }
  # get the offset
  offset <- function(y, w = 1)
  {
    if (!is.null(sigma)){
      RET <- FAM$sigma.linkfun(sigma)
    } else {
      eval(FAM$sigma.initial)
      RET <- FAM$sigma.linkfun(mean(sigma))
    }
    return(RET)
  }
  mboost:::Family(ngradient = ngradient, risk = risk, loss = loss,
         response = function(f) FAM$sigma.linkinv(f), offset = offset, 
         name = paste(FAM$family[2], "2nd parameter (sigma)"))
}

gamlss3parNu <- function(mu = NULL,  sigma = NULL, nu = NULL, fname="TF")
{
  FAM <- as.gamlss.family(fname)
  NAMEofFAMILY <- FAM$family
  dfun <- paste("d", fname, sep = "")
  pdf <- eval(parse(text = dfun))
  
  # get the loss
  loss <- function(y, f, w = 1, mu, sigma){
    - pdf(x = y, mu = mu, sigma = sigma,  nu = FAM$nu.linkinv(f),  log = TRUE) 
  }
  
  # compute the risk
  risk <- function(y, f, w = 1) {
    sum(w * loss(y = y, f = f, mu = mu, sigma = sigma))
  }
  # get the ngradient:  sigma is linkinv(f)
  # we need dl/deta = dl/dsigma*dsigma/deta
  ngradient <- function(y, f, w = 1){
    FAM$dldv(y = y, mu = mu, sigma = sigma, nu = FAM$nu.linkinv(f))*FAM$nu.dr(eta = f)
  }
  # get the offset
  offset <- function(y, w = 1)
  {
    if (!is.null(nu)){
      RET <- FAM$nu.linkfun(nu)
    } else {
      eval(FAM$nu.initial)
      RET <- FAM$nu.linkfun(mean(nu))
    }
    return(RET)
  }
  mboost:::Family(ngradient = ngradient, risk = risk, loss = loss,
         response = function(f) FAM$nu.linkinv(f), offset = offset, 
         name = paste(FAM$family[2], "3rd parameter (nu)"))
}



gamlss3parFam <- function(mu = NULL, sigma = NULL, nu = NULL,  fname = "TF")
{ 
  Families(mu = gamlss3parMu(mu = mu, sigma = sigma, nu = nu,  fname = fname), 
           sigma=gamlss3parSigma(mu = mu, sigma = sigma, nu = nu,  fname = fname),
           nu = gamlss3parNu(mu = mu, sigma = sigma, nu = nu,  fname = fname))
}

#--------------------------- 4 parameters 

gamlss4parMu <- function(mu = NULL,  sigma = NULL, nu = NULL, tau = NULL,  fname="BCP"){

  FAM <- as.gamlss.family(fname)
  NAMEofFAMILY <- FAM$family
  dfun <- paste("d", fname, sep = "")
  pdf <- eval(parse(text = dfun))
  
  # get the loss
  loss <- function(y, f, sigma, nu, tau,  w = 1){
    - pdf(x = y, mu = FAM$mu.linkinv(f), sigma = sigma, nu = nu, tau = tau,  log = TRUE) 
  }
  # compute the risk
  risk <- function(y, f, w = 1) {
    sum(w * loss(y = y, f = f, sigma = sigma, nu = nu, tau = tau))
  }
  # get the ngradient:  mu is linkinv(f)
  # we need dl/deta = dl/dmu*dmu/deta 
  ngradient <- function(y, f, w = 1){
    FAM$dldm(y = y, mu = FAM$mu.linkinv(f), sigma = sigma, nu = nu, tau = tau)*FAM$mu.dr(eta = f)
  }
  # get the offset -> we take the starting values of gamlss
  offset <- function(y, w = 1)
  {
    if (!is.null(mu)){
      RET <- FAM$mu.linkfun(mu)
    } else {
      eval(FAM$mu.initial)
      RET <- FAM$mu.linkfun(mean(mu))
    }
    return(RET)
  }
  mboost:::Family(ngradient = ngradient, risk = risk, loss = loss,
         response = function(f) FAM$mu.linkinv(f), offset=offset, 
         name = paste(FAM$family[2], "1st parameter (mu)"))
}


gamlss4parSigma <- function(mu = NULL,  sigma = NULL, nu = NULL, tau = NULL,  fname="BCPE")
{
  FAM <- as.gamlss.family(fname)
  NAMEofFAMILY <- FAM$family
  dfun <- paste("d", fname, sep = "")
  pdf <- eval(parse(text = dfun))
  
  # get the loss
  loss <- function(y, f, w = 1, mu, nu, tau){
    - pdf(x = y, mu = mu, sigma = FAM$sigma.linkinv(f), nu = nu, tau = tau,  log = TRUE) 
  }
  
  # compute the risk
  risk <- function(y, f, w = 1) {
    sum(w * loss(y = y, f = f, mu = mu, nu = nu, tau = tau))
  }
  # get the ngradient:  sigma is linkinv(f)
  # we need dl/deta = dl/dsigma*dsigma/deta
  ngradient <- function(y, f, w = 1){
    FAM$dldd(y = y, mu = mu, sigma = FAM$sigma.linkinv(f), nu = nu, tau = tau)*FAM$sigma.dr(eta = f)
  }
  # get the offset
  offset <- function(y, w = 1)
  {
    if (!is.null(sigma)){
      RET <- FAM$sigma.linkfun(sigma)
    } else {
      eval(FAM$sigma.initial)
      RET <- FAM$sigma.linkfun(mean(sigma))
    }
    return(RET)
  }
  mboost:::Family(ngradient = ngradient, risk = risk, loss = loss,
         response = function(f) FAM$sigma.linkinv(f), offset = offset, 
         name = paste(FAM$family[2], "2nd parameter (sigma)"))
}

gamlss4parNu <- function(mu = NULL,  sigma = NULL, nu = NULL, tau = NULL,  fname="BCPE")
{
  FAM <- as.gamlss.family(fname)
  NAMEofFAMILY <- FAM$family
  dfun <- paste("d", fname, sep = "")
  pdf <- eval(parse(text = dfun))
  
  # get the loss
  loss <- function(y, f, w = 1, mu, sigma, tau){
    - pdf(x = y, mu = mu, sigma = sigma,  nu = FAM$nu.linkinv(f), tau = tau,  log = TRUE) 
  }
  
  # compute the risk
  risk <- function(y, f, w = 1) {
    sum(w * loss(y = y, f = f, mu = mu, sigma = sigma, tau = tau))
  }
  # get the ngradient:  sigma is linkinv(f)
  # we need dl/deta = dl/dsigma*dsigma/deta
  ngradient <- function(y, f, w = 1){
    FAM$dldv(y = y, mu = mu, sigma = sigma, nu = FAM$nu.linkinv(f), tau  = tau)*FAM$nu.dr(eta = f)
  }
  # get the offset
  offset <- function(y, w = 1)
  {
    if (!is.null(nu)){
      RET <- FAM$nu.linkfun(nu)
    } else {
      eval(FAM$nu.initial)
      RET <- FAM$nu.linkfun(mean(nu))
    }
    return(RET)
  }
  mboost:::Family(ngradient = ngradient, risk = risk, loss = loss,
         response = function(f) FAM$nu.linkinv(f), offset = offset, 
         name = paste(FAM$family[2], "3rd parameter (nu)"))
}



gamlss4parTau <- function(mu = NULL,  sigma = NULL, nu = NULL, tau = NULL,  fname="BCPE")
{
  FAM <- as.gamlss.family(fname)
  NAMEofFAMILY <- FAM$family
  dfun <- paste("d", fname, sep = "")
  pdf <- eval(parse(text = dfun))
  
  # get the loss
  loss <- function(y, f, w = 1, mu, sigma, nu){
    - pdf(x = y, mu = mu, sigma = sigma,  nu =  nu , tau = FAM$tau.linkinv(f),  log = TRUE) 
  }
  
  # compute the risk
  risk <- function(y, f, w = 1) {
    sum(w * loss(y = y, f = f, mu = mu, sigma = sigma, nu = nu))
  }
  # get the ngradient
  ngradient <- function(y, f, w = 1){
    FAM$dldt(y = y, mu = mu, sigma = sigma, tau = FAM$tau.linkinv(f), nu  = nu)*FAM$tau.dr(eta = f)
  }
  # get the offset
  offset <- function(y, w = 1)
  {
    if (!is.null(tau)){
      RET <- FAM$tau.linkfun(tau)
    } else {
      eval(FAM$tau.initial)
      RET <- FAM$tau.linkfun(mean(tau))
    }
    return(RET)
  }
  mboost:::Family(ngradient = ngradient, risk = risk, loss = loss,
         response = function(f) FAM$tau.linkinv(f), offset = offset, 
         name = paste(FAM$family[2], "4th parameter (tau)"))
}


gamlss4parFam <- function(mu = NULL, sigma = NULL, nu = NULL, tau = NULL,  fname = "BCPE")
{ 

  Families(mu = gamlss4parMu(mu = mu, sigma = sigma, nu = nu, tau = tau,  fname = fname), 
           sigma=gamlss4parSigma(mu = mu, sigma = sigma, nu = nu, tau = tau,  fname = fname),
           nu = gamlss4parNu(mu = mu, sigma = sigma, nu = nu, tau = tau, fname = fname),
           tau = gamlss4parTau(mu = mu, sigma = sigma, nu = nu, tau = tau, fname = fname))
}


#------------- constructor
Families.gamlss <-function(fname = "NO",  mu = NULL, sigma = NULL, nu = NULL, tau = NULL)
  {
   # - require gamlss.dist
   if(! "gamlss.dist"  %in% library()$result) stop("Please install package 'gamlss.dist' for using gamlss families.")
   require(gamlss.dist, quietly = TRUE, warn.conflicts= FALSE)  
   
   family <- fname
   if (mode(family) != "character" && mode(family) != "name")
    fname <- as.character(substitute(family))
    
    npar <- eval(call(family))$nopar
    switch(npar,
         { # 1 parameter
           fun <- gamlss1parMu(mu = mu, fname = fname)
           warning("For boosting one-parametric families, please use the mboost package.")
           } ,
            { # 2 parameters
                    fun <- gamlss2parFam(mu = mu, sigma = sigma, fname = fname)
            },
            { # 3 parameters
                    fun <- gamlss3parFam(mu = mu, sigma = sigma, nu = nu , fname = fname)
                
            } ,
            { # 4 parameters
                    fun <- gamlss4parFam(mu = mu, sigma = sigma, nu = nu, tau = tau , fname = fname)
                
            } )
   fun
  }

