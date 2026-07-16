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
    },
    pcopula = function(u1, u2, theta){
      rho <- theta
      z1 <- qnorm(u1)
      z2 <- qnorm(u2)
      sapply(seq_along(u1), function(i)
      mvtnorm::pmvnorm(lower = c(-Inf, -Inf),
                              upper = c(z1[i], z2[i]),
                              corr = matrix(c(1, rho, rho, 1), nrow = 2)))
    },
    h = function(u1, u2, theta){
      rho <- theta
      z1 <- qnorm(u1); z2 <- qnorm(u2)
      pnorm((z1 - rho*z2) / sqrt(1 - rho^2))
    },
    response = function(f) tanh(f),
    dresponse = function(f) 1 - tanh(f)^2,
    response_inv = function(theta) atanh(theta),
    offset = 0,
    tau = function(theta) 2/pi * asin(theta),
    tau_inv = function(tau) sin(pi * tau/2),
    tau_range = c(-0.95, 0.95),
    name = "Gaussian copula: rho (tanh link)"
  )
}

.copula_clayton <- function(){
  list(
  pcopula = function(u1, u2, theta){
    (u1^(-theta) + u2^(-theta) - 1)^(-1/theta)
  },
  logdcopula = function(u1, u2, theta){
    log(theta + 1) - (theta + 1) * (log(u1) + log(u2)) - 
      (2 + 1/theta) * log(u1^(-theta) + u2^(-theta) - 1)
  },
  dlogdcopula_u1 = function(u1, u2, theta){
    term1 <- u1^(-theta) + u2^(-theta) -1 
    return(-(theta+1)/u1 + (2*theta+1) * u1^(-theta-1) / term1)
  },
  dlogdcopula_theta = function(u1, u2, theta){
    term1 <- u1^(-theta) + u2^(-theta) -1 
    return(1/(theta+1) - log(u1) - log(u2) + log(term1)/theta^2 + 
             (2 + 1/theta) * (u1^(-theta)*log(u1) + u2^(-theta)*log(u2)) / term1)
  },
  h = function(u1, u2, theta){
    u2^(-theta - 1) * (u1^(-theta) + u2^(-theta) - 1)^(-1 - 1/theta)
  },
  response = function(f) exp(f),
  dresponse = function(f) exp(f),
  response_inv = function(theta) log(theta),
  offset = 0,
  tau = function(theta) theta / (theta + 2),
  tau_inv = function(tau) 2 * tau / (1 - tau),
  tau_range = c(0.01, 0.95),
  name = "Clayton copula: theta (exp link)"
)
}

.copula_gumbel <- function(){
  list(
    pcopula = function(u1, u2, theta){
      x <- -log(u1); y <- -log(u2)
      term1 <- x^theta + y^theta
      exp(-term1^(1/theta))
    },
    logdcopula = function(u1, u2, theta){
      x <- -log(u1); y <- -log(u2)
      term1 <- x^theta + y^theta
      term2 <- exp(-term1^(1/theta))
      log(term2) - log(u1) - log(u2) + (2/theta - 2) * log(term1) +
        (theta - 1) * log(x * y) + log(1 + (theta - 1) * term1^(-1/theta))
    },
    dlogdcopula_u1 = function(u1, u2, theta){
      x <- -log(u1); y <- -log(u2)
      term1 <- x^theta + y^theta
      term2 <- term1^(1/theta)
      term3 <- 1 + (theta-1)/term2
      
      return((x^(theta-1)/term1 * (term2 + 2*theta - 2 + (theta-1)/
                                     (term2*term3)) - 1 - (theta-1)/x) / u1)
    },
    dlogdcopula_theta = function(u1, u2, theta){
      x <- -log(u1); y <- -log(u2)
      term1 <- x^theta + y^theta
      term12 <- x^theta * log(x) + y^theta * log(y)
      term2 <- term1^(1/theta)
      term3 <- 1 + (theta-1)/term2
      term4 <- log(term1)/theta^2 - term12/(theta*term1)
      
      return(term2*term4 - 2*log(term1)/theta^2 + (2/theta - 2)*term12/term1 +
               log(x*y) + (1 + (theta-1)*term4) / (term2 * term3))
    },
    h = function(u1, u2, theta){
      x <- -log(u1); y <- -log(u2)
      term1 <- x^theta + y^theta
      term2 <- exp(-term1^(1/theta))
      term2 * term1^(1/theta - 1) * y^(theta - 1) / u2
    },
    response = function(f) exp(f) + 1,
    dresponse = function(f) exp(f),
    response_inv = function(theta) log(theta - 1),
    offset = 0,
    tau = function(theta) 1 - 1/theta,
    tau_inv = function(tau) 1 / (1-tau),
    tau_range = c(0.01, 0.95),
    name = "Gumbel copula: theta (exp+1 link)"
  )
}

.copula_frank <- function(){
  list(
    pcopula = function(u1, u2, theta){
      term1 <- 1 - exp(-theta*u1)
      term2 <- 1 - exp(-theta*u2)
      -1/theta * log(1 - term1 * term2/(1 - exp(-theta)))
    },
    logdcopula = function(u1, u2, theta){
      term1 <- 1 - exp(-theta*u1)
      term2 <- 1 - exp(-theta*u2)
      D <- (1 - exp(-theta)) - term1 * term2
      log(theta * (1 - exp(-theta))) - theta*(u1+u2) - 2 * log(abs(D))
    },
    dlogdcopula_u1 = function(u1, u2, theta){
      a <- 1 - exp(-theta)
      t1 <- 1 - exp(-theta*u1)
      t2 <- 1 - exp(-theta*u2)
      D <- a - t1*t2
      theta * (-1 + 2*t2*exp(-theta*u1) / D)
    },
    dlogdcopula_theta = function(u1, u2, theta){
      a <- 1 - exp(-theta)
      t1 <- 1 - exp(-theta*u1)
      t2 <- 1 - exp(-theta*u2)
      D <- a - t1*t2
      dD <- exp(-theta) - u1*exp(-theta*u1)*t2 - t1*u2*exp(-theta*u2)
      1/theta + exp(-theta)/a - (u1+u2) - 2*dD/D
    },
    h = function(u1, u2, theta){
      (-exp(theta) * (exp(theta*u1) - 1)) / 
        (exp(theta*u1+theta*u2) - exp(theta*u1+theta) - exp(theta*u2 + theta) +
        exp(theta))
    },
    response = function(f) f,
    dresponse = function(f) 1,
    response_inv = function(theta) theta,
    offset = 1,
    tau = function(theta){      #no closed solution for tau-conversion
      sapply(theta, function(th){
        D1 <- integrate(function(t) ifelse(t == 0, 1, t/expm1(t)),
                        0, th)$value / th
        1 - 4/th * (1 - D1)
      })
    },
    tau_inv = function(tau){
      sapply(tau, function(tv){
        interval <- if (tv > 0) c(1e-6, 100) else c(-100, -1e-6)
        uniroot(function(th){
          D1 <- integrate(function(t) ifelse(t == 0, 1, t/expm1(t)),
                          0, th)$value / th
          1 - 4/th * (1-D1) - tv
        }, interval = interval, tol = 1e-9)$root
      })
    },
    tau_range = c(-0.95, 0.95),
    name = "Frank copula: theta (identity link)"
  )
}

# Input:  copula  - character with copula name
# Output: copula object 
get_copula <- function(copula){
  switch(copula,
         "gaussian" = .copula_gaussian(),
         "clayton" = .copula_clayton(),
         "gumbel" = .copula_gumbel(),
         "frank" = .copula_frank(),
         stop("Unknown copula: '", copula, "'"))
}
