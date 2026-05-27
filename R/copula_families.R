# Input:  marginal  - gamboostLSS-family-object
#         y         - response vector
#         params    - named list of current parameter values
#
#
#
# Output: cdf vector in (0,1) 
get_marginal_cdf <- function(marginal, y, params){
  # extract family name from family-object
  fam_name <- tolower(attr(marginal, "name"))
  
  # switch family to cdf
  switch(fam_name,
         "gaussian" = pnorm(y, mean = params$mu, sd = params$sigma),
         "gamma" = pgamma(y, shape = params$mu^2 / params$sigma^2,
                          rate = params$mu / params$sigma^2),
         "log-normal" = plnorm(y, meanlog = params$mu, sdlog = params$sigma),
         "beta" = pbeta(y, shape1 = params$mu * params$phi,
                        shape2 = (1 - params$mu) * params$phi),
         stop("unknown marginal family: '", fam_name, "'")
         )
}

# Input:  marginal    - families-object
#         y           - response vector
#         params      - list(mu = ..., sigma = ...)
#         param_name  - character: "mu" or "sigma" or "phi"
#
#
#
# Output: dF(y)/dparam - vector of same length as y 
dCDF_dparam <- function(marginal, y, params, param_name){
  # extract family name from family-object
  fam_name <- tolower(attr(marginal, "name"))
  
  # switch family to cdf derivative 
  # gaussian and log-normal analytic, Gamma, Beta numeric
  switch(fam_name, 
         "gaussian" = switch(param_name,
                             "mu" = -dnorm(y, mean = params$mu,
                                           sd = params$sigma),
                             "sigma" = -dnorm(y, mean = params$mu,
                                              sd = params$sigma) *
                               (y - params$mu) / params$sigma,
                             stop("Unknown parameter'", param_name, "' for
                                  gaussian")),
         "log-normal" = switch(param_name,
                               "mu" = -dnorm(log(y), mean = params$mu,
                                             sd = params$sigma),
                               "sigma" = -dnorm(log(y), mean = params$mu,
                                                sd = params$sigma) * 
                                 (log(y) - params$mu) / params$sigma,
                               stop("Unknown param '", param_name, "' for 
                                    log-normal")),
         "gamma" = ,
         "beta" = {
           eps <- 1e-6
           params_p <- params
           params_p[[param_name]] <- params[[param_name]] + eps
           params_m <- params
           params_m[[param_name]] <- params[[param_name]] - eps
           (get_marginal_cdf(marginal, y, params_p) -
               get_marginal_cdf(marginal, y, params_m)) / (2 * eps)
           },
         stop("Unknow family: '", fam_name, "'")
  )
}

# Input: 
#   sub_fam     - families-object 
#   f           - current predictor vector
#
# Output: dresponse(f)/df
get_dresponse <- function(sub_fam, f){
  eps <- 1e-7
  (sub_fam@response(f+eps) - sub_fam@response(f - eps)) / (2 * eps)
}




# Input: 
#   marginal1, marginal2  - family objects
#   copula                - character: "gaussian", "clayton", "gumbel", "frank"
#   theta                 - start value for copula-parameter
#
#
#
# Output: Family-object with named slots
#   mu1, sigma1 (Parameters of first marginal)
#   mu2, sigma2 (Parameters of second marginal)
#   theta (Copula-Parameter)
CopulaFamilies <- function(marginal1, marginal2, copula = "gaussian",
                           theta = 0, stabilization = c("none", "MAD", "L2")){
  
  stabilization <- check_stabilization()
  cop <- get_copula(copula)
  
  # Initialize the param start values
  pnames1 <- names(marginal1)
  pnames2 <- names(marginal2)
  slot_names1 <- paste0(pnames1, "1")
  slot_names2 <- paste0(pnames2, "2")
  
  # initialize closure variables
  for (nm in c(slot_names1, slot_names2)) assign(nm, NULL)
  theta_cop <- theta
  
  # Helper for params
  get_params <- function(){
    list(
      params1 = setNames(lapply(slot_names1, get), pnames1),
      params2 = setNames(lapply(slot_names2, get), pnames2)
    )
  }
  
  # Helper for marginal cdfs
  get_u <- function(y){
    p <- get_params()
    list(
      u1 = get_marginal_cdf(marginal1, y[,1], p$params1),
      u2 = get_marginal_cdf(marginal2, y[,2], p$params2)
    )
  }
  
  # Helper to validate copula response
  check_copula_response <- function(y){
    if(!is.matrix(y) || ncol(y) != 2)
      stop("Response must be a 2-column matrix for CopulaFamilies")
    return(y)
  }
  
  # Create the marginal family 
  marginal_family <- function(i, marginal, pnames, slot_names, side){
    local({
      pname <- pnames[i]
      sname <- slot_names[i]
      sub_fam <- marginal[[pname]]
      
      ngradient_fn <- function(y, f, w = 1){
        p <- get_params()
        u1 <- get_marginal_cdf(marginal1, y[,1], p$params1)
        u2 <- get_marginal_cdf(marginal2, y[,2], p$params2)
        y_m <- y[, side]
        pm <- if (side == 1) p$params1 else p$params2
        marg <- if (side == 1) marginal1 else marginal2
        u_act <- if (side == 1) u1 else u2
        u_oth <- if (side == 1) u2 else u1
        
        for (other_param in pnames[pnames != pname])
          assign(other_param, pm[[other_param]], environment(sub_fam@ngradient))
        ngr_marginal <- sub_fam@ngradient(y_m, f, w = 1)
        
        # copula correction
        cop_corr <- cop$dlogdcopula_u1(u_act, u_oth, theta_cop) * 
          dCDF_dparam(marg, y_m, pm, pname) * get_dresponse(sub_fam, f)
        
        stabilize_ngradient(ngr_marginal + cop_corr, w = w, stabilization)
      }
      
      loss_function <- function(y, f){
        u <- get_u(y) 
        return(- cop$logdcopula(u$u1, u$u2, theta_cop))
      }
      
      offset_fn <- function(y, w){
        return(sub_fam@offset(y[, side], w))
      }
      
      # build family object
      Family(
        ngradient = ngradient_fn,
        loss = loss_function,
        risk = function(y, f, w = 1) sum(w * loss_function(y, f)),
        response = sub_fam@response,
        offset = offset_fn, 
        check_y = check_copula_response,
        name = paste0("CopulaFamilies: ", sname)
      )
    })
  }
  
  
  
}











