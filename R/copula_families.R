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
                           theta = 0){
  cop <- get_copula(copula)
  
  marginal_type_1 <- get_marginal_type(marginal1)
  marginal_type_2 <- get_marginal_type(marginal2)
  
}