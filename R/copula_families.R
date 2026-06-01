# Input:  marginal  - families-object
# Output: character "continuous", "discrete" or "binary"
get_marginal_type <- function(marginal){
  fam_name <- tolower(attr(marginal, "name"))
  
  continuous <- c("gaussian", "gamma", "log-normal", "beta", "student t",
                  "log-log", "weibull")
  discrete <- c("zip", "negative binomial", "zinbi")
  binary <- c(0)
  
  if (fam_name %in% continuous) return("continuous")
  if (fam_name %in% discrete) return("discrete")
  if (fam_name %in% binary) return("binary")
  
  stop("Unknown marginal type for family: '", fam_name, "'")
}

# Input:  marginal1, marginal2   - families-object
# Output: character, e.g. "continuous_continuous", "binary_continuous",...
get_marginal_case <- function(marginal1, marginal2){
  type1 <- get_marginal_type(marginal1)
  type2 <- get_marginal_type(marginal2)
  
  paste0(type1, "_", type2)
}

# Input:    marginal  - families-object
#           y         - response vector
#           params    - named list of current parameter values
# Output:   F(y|params) - vector in (0,1)
get_marginal_cdf <- function(marginal, y, params){
  fam_name <- tolower(attr(marginal, "name"))
  
  switch(fam_name,
         "gaussian" = pnorm(y, mean = params$mu, sd = params$sigma),
         "gamma" = pgamma(y, shape = params$mu^2 / params$sigma^2,
                          rate = params$mu / params$sigma^2),
         "log-normal" = plnorm(y, meanlog = params$mu, sdlog = params$sigma),
         "beta" = pbeta(y, shape1 = params$mu * params$phi,
                        shape2 = (1-params$mu) * params$phi),
         stop("Unknown family without cdf: '", fam_name, "'"))
}

# Input:    marginal  - families-object
#           y         - response vector
#           params    - named list of current parameter values
# Output:   log f(y|params) - vector
get_marginal_logpdf <- function(marginal, y, params){
  fam_name <- tolower(attr(marginal, "name"))
  
  switch(fam_name,
         "gaussian" = dnorm(y, mean = params$mu, sd = params$sigma, log = T),
         "gamma" = dgamma(y, shape = params$mu^2 / params$sigma^2,
                          rate = params$mu / params$sigma^2, log = T),
         "log-normal" = dlnorm(y, meanlog = params$mu, sdlog = params$sigma,
                               log = T),
         "beta" = dbeta(y, shape1 = params$mu * params$phi,
                        shape2 = (1-params$mu) * params$phi, log = T),
         stop("Unknown family without pdf: '", fam_name, "'"))
}

# Input:    y                       - n x 2 matrix of observations
#           marginal1, marginal2    - families-objects
#           params1, params2        - named lists of parameter values
#           copula                  - function(u1, u2, theta) 
#           theta                   - copula dependence parameter
# Output:   vector of length n with individual log-likelihood contributions
continuous_continuous_loglik <- function(y, marginal1, marginal2, params1,
                                         params2, copula, theta){
  
  cop <- get_copula(copula) # get_copula will be implemented tomorrow in 
                            # new script along with storing the copula 
                            # properties for each copula 
  u1 <- get_marginal_cdf(marginal1, y[,1], params1)
  u2 <- get_marginal_cdf(marginal2, y[,2], params2)
  
  cop$logdcopula(u1, u2, theta) + 
    get_marginal_logpdf(marginal1, y[,1], params1) +
    get_marginal_logpdf(marginal2, y[,2], params2)
}




