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
         "gamma" = pgamma(y, shape = params$sigma,
                          rate = params$sigma / params$mu),
         "log-normal" = plnorm(y, meanlog = params$mu, sdlog = params$sigma),
         "beta" = pbeta(y, shape1 = params$mu * params$phi,
                        shape2 = (1-params$mu) * params$phi),
         "student t" = pt((y - params$mu) / params$sigma, df = params$df),
         "weibull" = pweibull(y, shape = 1/params$sigma, 
                              scale = exp(params$mu)),
         "log-log" = plogis((log(y) - params$mu) / params$sigma),
         "negative binomial" = pnbinom(y, size = params$sigma, mu = params$mu),
         "zip" = params$sigma + (1-params$sigma) * ppois(y, lambda = params$mu),
         "zinbi" = params$nu + (1-params$nu) * 
           pnbinom(y, size = 1/params$sigma, mu = params$mu),
         stop("Unknown family without cdf: '", fam_name, "'"))
}

# Input:    marginal  - families-object
#           y         - response vector
#           params    - named list of current parameter values
# Output:   log f(y|params) - vector
get_marginal_logpdf <- function(marginal, y, params){
  fam_name <- tolower(attr(marginal, "name"))
  
  switch(fam_name,
         "gaussian" = dnorm(y, mean = params$mu, sd = params$sigma, log = TRUE),
         "gamma" = dgamma(y, shape = params$sigma,
                          rate = params$sigma / params$mu, log = TRUE),
         "log-normal" = dlnorm(y, meanlog = params$mu, sdlog = params$sigma,
                               log = TRUE),
         "beta" = dbeta(y, shape1 = params$mu * params$phi,
                        shape2 = (1-params$mu) * params$phi, log = TRUE),
         "student t" = -log(params$sigma) + dt((y - params$mu) / params$sigma,
                                               df = params$df, log = TRUE),
         "weibull" = dweibull(y, shape = 1/params$sigma, 
                              scale = exp(params$mu), log = TRUE),
         "log-log" = dlogis((log(y) - params$mu) / params$sigma, log = TRUE) - 
           log(params$sigma) - log(y),
         "negative binomial" = dnbinom(y, size = params$sigma, mu = params$mu,
                                       log = TRUE),
         "zip" = ifelse(y == 0, 
                        log(params$sigma + (1-params$sigma) * 
                          dpois(y, lambda = params$mu)),
                        log((1-params$sigma) * 
                          dpois(y, lambda = params$mu))),
         "zinbi" = ifelse(y == 0, log(params$nu + (1-params$nu) * 
                            dnbinom(y, size = 1/params$sigma, mu = params$mu)),
                          log((1-params$nu) * 
                            dnbinom(y, size = 1/params$sigma, mu = params$mu))),
         stop("Unknown family without pdf: '", fam_name, "'"))
}

# Input:    y                       - n x 2 matrix of observations
#           marginal1, marginal2    - families-objects
#           params1, params2        - named lists of parameter values
#           copula                  - character, e.g. "gaussian"
#           theta                   - copula dependence parameter
# Output:   vector of length n with individual log-likelihood contributions
continuous_continuous_loglik <- function(y, marginal1, marginal2, params1,
                                         params2, copula, theta){
  
  cop <- get_copula(copula) 
  u1 <- get_marginal_cdf(marginal1, y[,1], params1)
  u2 <- get_marginal_cdf(marginal2, y[,2], params2)
  
  cop$logdcopula(u1, u2, theta) + 
    get_marginal_logpdf(marginal1, y[,1], params1) +
    get_marginal_logpdf(marginal2, y[,2], params2)
}

# Input:    y                       - n x 2 matrix of observations
#           marginal1, marginal2    - families-objects
#           params1, params2        - named lists of parameter values
#           copula                  - character, e.g. "gaussian" 
#           theta                   - copula dependence parameter
# Output:   vector of length n with individual log-likelihood contributions
discrete_discrete_loglik <- function(y, marginal1, marginal2, params1,
                                     params2, copula, theta){
  
  cop <- get_copula(copula) 
  u1 <- get_marginal_cdf(marginal1, y[,1], params1)
  u2 <- get_marginal_cdf(marginal2, y[,2], params2)
  u1_minus <- get_marginal_cdf(marginal1, y[,1] - 1, params1)
  u2_minus <- get_marginal_cdf(marginal2, y[,2] - 1, params2)
  
  p <- cop$pcopula(u1, u2, theta) - cop$pcopula(u1_minus, u2, theta) -
    cop$pcopula(u1, u2_minus, theta) + cop$pcopula(u1_minus, u2_minus, theta)
  
  p <- pmax(p, 1e-10)
  return(log(p))
}

# Input:    marginal1, marginal2  - gamboostLSS families-objects
#           copula                - character e.g. "gaussian"
#           theta                 - starting value for copula parameter
#           stabilization         - character e.g. "MAD" 
# Output:   families-object
CopulaFamilies <- function(marginal1, marginal2, copula = "gaussian", theta = 1,
                           stabilization = c("none", "MAD", "L2")){
  
  # Identify the correct copula first 
  cop <- get_copula(copula)
  
  # Identify the case and switch 
  marginal_case <- get_marginal_case(marginal1, marginal2)
  loglik <- switch(marginal_case,
         "continuous_continuous" = continuous_continuous_loglik,
         "discrete_discrete" = discrete_discrete_loglik)
  
  # Use stabilization-method
  stabilization <- check_stabilization(stabilization)
  
  # Extract the names of the marginals
  params1 <- names(marginal1)
  params2 <- names(marginal2)
  
  # Extract the quantile functions
  qfun1 <- attr(marginal1, "qfun")
  qfun2 <- attr(marginal2, "qfun")
  
  # closure variables
  current_params1 <- setNames(vector("list", length(params1)), params1)
  current_params2 <- setNames(vector("list", length(params2)), params2)
  
  
  # marginal 1
  sub_families1 <- list()
  for (p in params1){
    sub_families1[[p]] <- local({
      param_name <- p
      
      # loss 
      loss_p <- function(y, f){
        current_params1[[param_name]] <- marginal1[[param_name]]@response(f)
        return(-loglik(y, marginal1, marginal2, current_params1, 
                              current_params2, copula, theta))
      }
      
      # ngradient
      ngradient_p <- function(y, f, w=1){NULL}
      
      # offset
      offset_p <- function(y, w){
        RET <- marginal1[[param_name]]@offset(y[,1], w)
        current_params1[[param_name]] <<- marginal1[[param_name]]@response(RET)
        return(RET)
      }
      
      Family(loss = loss_p,
             ngradient = ngradient_p,
             offset = offset_p,
             response = marginal1[[param_name]]@response,
             name = marginal1[[param_name]]@name)
    })
  }
  
  # marginal 2
  sub_families2 <- list()
  for (p in params2){
    sub_families2[[p]] <- local({
      param_name <- p
      
      # loss
      loss_p <- function(y, f){
        current_params2[[param_name]] <- marginal2[[param_name]]@response(f)
        return(-loglik(y, marginal1, marginal2, current_params1, 
                       current_params2, copula, theta))
      }
      
      # ngradient
      ngradient_p <- function(y, f, w=1){NULL}
      
      # offset 
      offset_p <- function(y, w){
        RET <- marginal2[[param_name]]@offset(y[,2], w)
        current_params2[[param_name]] <<- marginal2[[param_name]]@response(RET)
        return(RET)
      }
      
      Family(loss = loss_p,
             ngradient = ngradient_p,
             offset = offset_p,
             response = marginal2[[param_name]]@response,
             name = marginal2[[param_name]]@name)
    })
  }
  
  # copula
  # loss
  loss_theta <- function(y, f){
    theta_current <- f # no response for dependence param implemented yet
    return(-loglik(y, marginal1, marginal2, current_params1, 
                   current_params2, copula, theta_current))
  }
  
  # ngradient
  ngradient_theta <- function(y, f, w=1){
    NULL
    }
  
  # offset
  offset_theta <- function(y, w){
    NULL
  }
  # also not possible to implement yet, since we need response-function for 
  # chosen copula
  
  theta_family <- Family(loss = loss_theta,
         ngradient = ngradient_theta,
         offset = offset_theta,
         response = cop$response,
         name = cop$name
         )
  
  # combine the individual family-objects into one for joint family
  do.call(Families, c(sub_families1, sub_families2, 
                      list(theta = theta_family),
                      list(qfun = NULL, 
                           name = paste0("Copula Families: ", copula))
                      ))
}











