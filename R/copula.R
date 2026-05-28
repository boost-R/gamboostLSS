.copula_gaussian <- function(){
  list(
    logdcopula = function(u1, u2, theta){
      rho <- theta
      z1 <- qnorm(u1)
      z2 <- qnorm(u2)
      
      return(-0.5*log(1 - rho^2) + (2*rho*z1*z2 - rho^2*(z1^2 + z2^2)) / 
               (2 * (1-rho^2)))
    },
    dlogdcopula_u1 = function(u1, u2, theta){
      rho <- theta
      z1 <- qnorm(u1)
      z2 <- qnorm(u2)
      
      return((rho*z2 - rho^2*z1) / ((1-rho^2) * dnorm(z1)))
    },
    dlogdcopula_theta = function(u1, u2, theta){
      rho <- theta
      z1 <- qnorm(u1)
      z2 <- qnorm(u2)
      
      return(rho / (1 - rho^2) + (z1^2 + z2^2 - 2*rho*z1*z2) * rho /
               (1-rho^2)^2 - z1*z2 / (1 - rho^2))
    }
  )
}

.copula_clayton <- function(){
  list(
    logdcopula = function(u1, u2, theta){
      return(log(1 + theta) + (-1 - theta)*(log(u1) + log(u2)) + 
               (-1/theta - 2) * log(u1^(-theta) + u2^(-theta) - 1))
    },
    dlogdcopula_u1 = function(u1, u2, theta){
      return((-1 - theta)/u1 + (-1/theta - 2) * (-theta * u1^(-theta-1)) / 
               (u1^(-theta) + u2^(-theta) -1))
    },
    dlogdcopula_theta = function(u1, u2, theta){
      term1 <- u1^(-theta) + u2^(-theta) - 1
      
      return(1 / (1 + theta) + 
               (log(u1) + log(u2)) +
               (1/theta^2) * log(term1) + 
               (-1/theta - 2) * (-(u1^(-theta) * log(u1) + u2^(-theta) * 
                                     log(u2))) /term1)
    }
  )
}

.copula_gumbel <- function(){
  list(
    logdcopula = function(u1, u2, theta){
      t1 <- (-log(u1))^theta
      t2 <- (-log(u2))^theta
      t <- t1 + t2
      C <- exp(-t^(1/theta))
      
      return(log(C) + (1/theta - 2) * log(t) + 
               log(t1*t2) - log(u1*u2) - log((-log(u1))*(-log(u2))) +
               log(t^(1/theta) + theta - 1)
             )
    },
    dlogdcopula_u1 = function(u1, u2, theta){
      a1 <- -log(u1); a2 <- -log(u2)
      t1 <- a1^theta; t2 <- a2^theta
      t <- t1 + t2 
      s <- t^(1/theta)
      
      dt_du1 <- -theta * a1^(theta -1) / u1
      
      term_C <- a1^(theta-1) * t^(1/theta - 1) / u1
      term_t <- (1/theta - 2) * dt_du1 / t
      term_a1 <- -(theta - 1) / (a1 * u1)
      term_logu <- -1/u1
      ds_du1 <- (1/theta) * t^(1/theta - 1) * dt_du1
      term_end <- ds_du1 / (s + theta - 1)
      
      return(term_C + term_t + term_a1 + term_logu + term_end)
      
    },
    dlogdcopula_theta = function(u1, u2, theta){
      a1 <- -log(u1); a2 <- -log(u2)
      t1 <- a1^theta; t2 <- a2^theta
      t <- t1 + t2
      s <- t^(1/theta)
      
      dt_theta <- t1*log(a1) + t2*log(a2)
      ds_theta <- s * (-log(t)/theta^2 + dt_theta/(theta*t))
      
      term_C <- -ds_theta
      term_t <- -1/theta^2 * log(t) + (1/theta - 2)*dt_theta/t
      term_a <- log(a1) + log(a2)
      term_end <- (ds_theta + 1) / (s + theta -1)
      
      return(term_C + term_t + term_a + term_end)
    }
  )
}

.copula_frank <- function(){
  list(
    logdcopula = function(u1, u2, theta){
      a <- expm1(-theta * u1)
      b <- expm1(-theta * u2)
      D <- expm1(-theta) + a * b
      
      return(log(abs(theta)) + log(abs(expm1(-theta))) - theta*(u1 + u2) -
               2*log(abs(D)))
    },
    dlogdcopula_u1 = function(u1, u2, theta){
      a <- expm1(-theta * u1)
      b <- expm1(-theta * u2)
      D <- expm1(-theta) + a * b
      
      return(theta * (2*b*(1+a)/D -1))
    },
    dlogdcopula_theta = function(u1, u2, theta){
      eps <- 1e-6
      cop <- .copula_frank()
      
      return((cop$logdcopula(u1, u2, theta+eps) -
               cop$logdcopula(u1, u2, theta - eps)) / (2*eps))
    }
  )
}

get_copula <- function(copula){
  switch(copula, 
         gaussian = .copula_gaussian(),
         clayton = .copula_clayton(),
         gumbel = .copula_gumbel(),
         frank = .copula_frank(),
         stop("Unknown copula: '", copula, "'. Choose from: gaussian,
              clayton, gumbel or frank.")
         )
}