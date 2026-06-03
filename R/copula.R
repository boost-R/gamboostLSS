.copula_gaussian <- function(){
  list(
    logdcopula = function(u1, u2, theta){
      rho <- theta
      z1 <- qnorm(u1)
      z2 <- qnorm(u2)
      return(-0.5 * log(1 - rho^2) + 
               (2*rho*z1*z2 - rho^2*(z1^2 + z2^2)) / (2*(1-rho^2))
      )
    },
    dlogdcopula_u1 = function(u1, u2, theta){
      rho <- theta
      z1 <- qnorm(u1)
      z2 <- qnorm(u2)
      return((rho*z2 - rho^2*z1) / ((1 - rho^2) * dnorm(z1)))
    },
    dlogdcopula_theta = function(u1, u2, theta){
      rho <- theta
      z1 <- qnorm(u1)
      z2 <- qnorm(u2)
      return(
        rho / (1 - rho^2) + (z1*z2*(1+rho^2) - rho*(z1^2+z2^2)) / (1 - rho^2)^2
      )
    }
  )
}

# Input:  copula  - character with copula name
# Output: copula object 
get_copula <- function(copula){
  switch(copula,
         "gaussian" = .copula_gaussian(),
         stop("Unknown copula: '", copula, "'"))
}
