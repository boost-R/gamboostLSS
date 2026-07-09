# Input:  marginal  - families-object
# Output: character "continuous", "discrete" or "binary"
get_marginal_type <- function(marginal){
  fam_name <- tolower(attr(marginal, "name"))
  
  continuous <- c("gaussian", "gamma", "log-normal", "beta", "student t",
                  "log-log", "weibull")
  discrete <- c("zip", "negative binomial", "zinbi")
  binary <- c("bernoulli (probit)", "bernoulli (logit)", "bernoulli (cloglog)")
  
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
         "bernoulli (probit)" = ,
         "bernoulli (logit)" = ,
         "bernoulli (cloglog)" = ifelse(y<0, 0, ifelse(y<1, 1 - params$mu, 1)),
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

# Input:    y                       - n x 2 matrix of responses
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

# Input:    y                       - n x 2 matrix of responses
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

# Input:    y                       - n x 2 matrix of responses
#           marginal1, marginal2    - families-objects
#           params1, params2        - named lists of parameter values
#           copula                  - character, e.g. "gaussian" 
#           theta                   - copula dependence parameter
# Output:   vector of length n with individual log-likelihood contributions
binary_continuous_loglik <- function(y, marginal1, marginal2, params1, 
                                     params2, copula, theta){
  cop <- get_copula(copula)
  u1 <- get_marginal_cdf(marginal1, y[,1], params1)
  u1_minus <- get_marginal_cdf(marginal1, y[,1]-1, params1)
  u2 <- get_marginal_cdf(marginal2, y[,2], params2)
  
  d <- cop$h(u1, u2, theta) - cop$h(u1_minus, u2, theta)
  d <- pmax(d, 1e-10)
  return(log(d) + get_marginal_logpdf(marginal2, y[,2], params2))
  
}

# Input:    y                       - n x 2 matrix of responses
#           marginal1, marginal2    - families-objects
#           params1, params2        - named lists of parameter values
#           copula                  - character, e.g. "gaussian" 
#           theta                   - copula dependence parameter
# Output:   vector of length n with individual log-likelihood contributions
continuous_binary_loglik <- function(y, marginal1, marginal2, params1, 
                                     params2, copula, theta){
  cop <- get_copula(copula)
  u1 <- get_marginal_cdf(marginal1, y[,1], params1)
  u2 <- get_marginal_cdf(marginal2, y[,2], params2)
  u2_minus <- get_marginal_cdf(marginal2, y[,2]-1, params2)
  
  d <- cop$h(u2, u1, theta) - cop$h(u2_minus, u1, theta)  
  d <- pmax(d, 1e-10)
  return(log(d) + get_marginal_logpdf(marginal1, y[,1], params1))
  
}

# Input:    marginal              - gamboostLSS families-object
#           y                     - response vector
#           params                - named list of parameter values
#           param_name            - character e.g. "mu"
# Output    
get_marginal_dcdf <- function(marginal, y, params, param_name){
  fam_name <- tolower(attr(marginal, "name"))
  switch(fam_name,
         "gaussian" = {
           if (param_name == "mu")
             -dnorm(y, mean = params$mu, sd = params$sigma)
           else if (param_name == "sigma")
             -dnorm(y, mean = params$mu, sd = params$sigma) * (y - params$mu) /
             params$sigma
         },
         sapply(y, function(yi){
           numDeriv::grad(function(p){
             params_p <- params
             params_p[[param_name]] <- p
             get_marginal_cdf(marginal, yi, params_p)
           }, params[[param_name]])
         })
  )
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
         "discrete_discrete" = discrete_discrete_loglik,
         "binary_binary" = discrete_discrete_loglik,
         "binary_continuous" = binary_continuous_loglik,
         "continuous_binary" = continuous_binary_loglik)
  
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
  current_theta <- theta
  
  update_current_params <- function(env){
    for (p in params1){
      nm <- paste0(p, "1")
      if (exists(nm, envir = env, inherits = FALSE))
        current_params1[[p]] <<- get(nm, envir = env)
    }
    for (p in params2){
      nm <- paste0(p, "2")
      if (exists(nm, envir = env, inherits = FALSE))
        current_params2[[p]] <<- get(nm, envir = env)
    }
    if (exists("theta", envir = env, inherits = FALSE))
      current_theta <<- get("theta", envir = env)
  }
  
  # marginal 1
  sub_families1 <- list()
  for (p in params1){
    sub_families1[[p]] <- local({
      param_name <- p
      
      # loss 
      loss_p <- function(y, f){
        current_params1[[param_name]] <- marginal1[[param_name]]@response(f)
        return(-loglik(y, marginal1, marginal2, current_params1, 
                              current_params2, copula, current_theta))
      }
      
      # ngradient
      ngradient_p <- function(y, f, w=1){
        update_current_params(parent.env(environment()))
        params1_f <- current_params1
        params1_f[[param_name]] <- marginal1[[param_name]]@response(f)
        
        eps <- 1e-7
        pv <- params1_f[[param_name]]
        p_plus <- params1_f; p_plus[[param_name]] <- pv + eps
        p_minus <- params1_f; p_minus[[param_name]] <- pv - eps
        dlogf1_dparam <- (get_marginal_logpdf(marginal1, y[,1], p_plus) -
                            get_marginal_logpdf(marginal1, y[,1], p_minus)) / 
                            (2*eps)
        dresponse_f <- sapply(f, function(fi) 
                              numDeriv::grad(marginal1[[param_name]]@response,
                                             fi))
        marg_ngradient <- dlogf1_dparam * dresponse_f
        
        if (marginal_case != "continuous_continuous" ||
            is.null(cop$dlogdcopula_u1))
          return(marg_ngradient)
        
        u1 <- get_marginal_cdf(marginal1, y[,1], params1_f)
        u2 <- get_marginal_cdf(marginal2, y[,2], current_params2)
        dcdf <- get_marginal_dcdf(marginal1, y[,1], params1_f, param_name)
        
        return(cop$dlogdcopula_u1(u1, u2, current_theta) * dcdf * dresponse_f + 
                 marg_ngradient)
      }
      
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
  names(sub_families1) <- paste0(params1, "1")
  
  # marginal 2
  sub_families2 <- list()
  for (p in params2){
    sub_families2[[p]] <- local({
      param_name <- p
      
      # loss
      loss_p <- function(y, f){
        current_params2[[param_name]] <- marginal2[[param_name]]@response(f)
        return(-loglik(y, marginal1, marginal2, current_params1, 
                       current_params2, copula, current_theta))
      }
      
      # ngradient
      ngradient_p <- function(y, f, w=1){
        update_current_params(parent.env(environment()))
        params2_f <- current_params2
        params2_f[[param_name]] <- marginal2[[param_name]]@response(f)
        
        eps <- 1e-7
        pv <- params2_f[[param_name]]
        p_plus <- params2_f; p_plus[[param_name]] <- pv + eps
        p_minus <- params2_f; p_minus[[param_name]] <- pv - eps
        dlogf2_dparam <- (get_marginal_logpdf(marginal2, y[,2], p_plus) -
                            get_marginal_logpdf(marginal2, y[,2], p_minus)) / 
                            (2*eps)
        dresponse_f <- sapply(f, function(fi) 
          numDeriv::grad(marginal2[[param_name]]@response,
                         fi))
        marg_ngradient <- dlogf2_dparam * dresponse_f
        
        if (marginal_case != "continuous_continuous" ||
            is.null(cop$dlogdcopula_u1))
          return(marg_ngradient)
        
        u1 <- get_marginal_cdf(marginal1, y[,1], current_params1)
        u2 <- get_marginal_cdf(marginal2, y[,2], params2_f)
        dcdf <- get_marginal_dcdf(marginal2, y[,2], params2_f, param_name)
        
        return(cop$dlogdcopula_u1(u1, u2, current_theta) * dcdf * dresponse_f + 
                 marg_ngradient)
      }
      
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
  
  names(sub_families2) <- paste0(params2, "2")
  # copula
  # loss
  loss_theta <- function(y, f){
    theta_current <- cop$response(f) 
    return(-loglik(y, marginal1, marginal2, current_params1, 
                   current_params2, copula, theta_current))
  }
  
  # ngradient
  ngradient_theta <- function(y, f, w=1){
    theta_current <- cop$response(f)
    u1 <- get_marginal_cdf(marginal1, y[,1], current_params1)
    u2 <- get_marginal_cdf(marginal2, y[,2], current_params2)
    outer_deriv <- cop$dlogdcopula_theta(u1, u2, theta_current)
    inner_deriv <- cop$dresponse(f)
    return(outer_deriv * inner_deriv)
    }
  
  # offset
  offset_theta <- function(y, w){
    return(0)
  }

  
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

