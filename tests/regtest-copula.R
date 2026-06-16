require("gamboostLSS")

### check continuous_continuous_loglik() accuracy at theta = 0
theta <- 0
y <- cbind(c(-2, -1, 0, 1, 2), c(-2, -1, 0, 1, 2))
stopifnot(max(abs(gamboostLSS:::continuous_continuous_loglik(y, 
                                                             GaussianLSS(), 
                                                             GaussianLSS(), 
                                                             list(mu=0, sigma=1), 
                                                             list(mu=0, sigma=1), 
                                                             "gaussian", 0) -
                    (gamboostLSS:::get_marginal_logpdf(GaussianLSS(),
                                                       y[,1],
                                                       list(mu=0, sigma=1)) +
                       gamboostLSS:::get_marginal_logpdf(GaussianLSS(),
                                                         y[,2],
                                                         list(mu=0, sigma=1)
                                                         )))) < 1e-6)

### check discrete_discrete_loglik() accuracy at theta = 0
theta <- 0
y <- cbind(c(0, 1, 2, 5), c(0, 2, 3, 1))

stopifnot(max(abs(gamboostLSS:::discrete_discrete_loglik(y, 
                                                             NBinomialLSS(), 
                                                             NBinomialLSS(), 
                                                             list(mu=2, sigma=1), 
                                                             list(mu=2, sigma=1), 
                                                             "gaussian", theta) -
                    (gamboostLSS:::get_marginal_logpdf(NBinomialLSS(),
                                                       y[,1],
                                                       list(mu=2, sigma=1)) +
                       gamboostLSS:::get_marginal_logpdf(NBinomialLSS(),
                                                         y[,2],
                                                         list(mu=2, sigma=1)
                       )))) < 1e-3) # lower accuracy threshold due to accuracy 
                                    # of pmvnorm

### check discrete_discrete_loglik() accuracy for large theta and extreme 
### combination of response
theta <- 0.95
y <- cbind(0, 20)
stopifnot(is.finite(gamboostLSS:::discrete_discrete_loglik(y, 
                                                           NBinomialLSS(), 
                                                           NBinomialLSS(), 
                                                           list(mu=2, sigma=1), 
                                                           list(mu=2, sigma=1), 
                                                           "gaussian", theta)))

### numerical gradient check
cop <- gamboostLSS:::get_copula("gaussian")
test_points <- list(
  # standard case - all params non-problematic
  list(u1 = 0.3, u2 = 0.6, theta = 0.4), 
  # negative correlation
  list(u1 = 0.7, u2 = 0.2, theta = -0.5),
  # strong dependence 
  list(u1 = 0.5, u2 = 0.5, theta = 0.85),
  # no correlation
  list(u1 = 0.1, u2 = 0.9, theta = 0)
)
eps <- 1e-6

# check analytical against numerical gradient of gaussian-copula (u1)
for (p in test_points){
  num_grad_u1 <- (cop$logdcopula(p$u1 + eps, p$u2, p$theta) - 
                    cop$logdcopula(p$u1 - eps, p$u2, p$theta)) / (2*eps)
  ana_grad_u1 <- cop$dlogdcopula_u1(p$u1, p$u2, p$theta)
  
  stopifnot(abs(num_grad_u1 - ana_grad_u1) < 1e-5)
}

# check analytical against numerical gradient of gaussian-copula (theta)
for (p in test_points){
  num_grad_theta <- (cop$logdcopula(p$u1, p$u2, p$theta + eps) -
                       cop$logdcopula(p$u1, p$u2, p$theta - eps)) / (2*eps)
  ana_grad_theta <- cop$dlogdcopula_theta(p$u1, p$u2, p$theta)
  
  stopifnot(abs(num_grad_theta - ana_grad_theta) < 1e-5)
}

### check numerical gradient of continuous_continuous_loglik
set.seed(42)
y <- cbind(rlnorm(10), rbeta(10, 5, 1))
marginal1 <- LogNormalLSS()
marginal2 <- BetaLSS()
params1 <- list(mu = 0.25, sigma = 0.95)
params2 <- list(mu = 0.8, phi = 5)
copula <- "gaussian"
theta <- 0.5

eps <- 1e-6
u1 <- gamboostLSS:::get_marginal_cdf(marginal1, y[,1], params1)
u2 <- gamboostLSS:::get_marginal_cdf(marginal2, y[,2], params2)
cop <- gamboostLSS:::get_copula(copula)

num_grad_loglik <- (gamboostLSS:::continuous_continuous_loglik(y, 
                                                               marginal1,
                                                               marginal2,
                                                               params1,
                                                               params2,
                                                               copula,
                                                               theta + eps) - 
                      gamboostLSS:::continuous_continuous_loglik(y, 
                                                                 marginal1,
                                                                 marginal2,
                                                                 params1,
                                                                 params2,
                                                                 copula,
                                                                 theta - eps)) /
  (2*eps)
ana_grad_loglik <- cop$dlogdcopula_theta(u1, u2, theta)

stopifnot(max(abs(num_grad_loglik - ana_grad_loglik)) < 1e-5)

### reference value test for discrete_discrete_loglik 
y <- cbind(1, 2)
theta <- 0.5
stopifnot(abs(gamboostLSS:::discrete_discrete_loglik(y, 
                                                     NBinomialLSS(),
                                                     NBinomialLSS(),
                                                     list(mu=2, sigma=1),
                                                     list(mu=2, sigma=1),
                                                     "gaussian", theta) - 
                (-3.328349)) < 1e-3)  # reference computed manually via 
                                      # mvtnorm::pmvnorm 

### check h against numeric derivative of pcopula 
set.seed(42)
y <- cbind(rnorm(10), rnorm(10))
eps <- 1e-3
copula <- "gaussian"
marginal1 <- GaussianLSS()
marginal2 <- GaussianLSS()
params1 <- list(mu = 0, sigma = 1)
params2 <- list(mu = 0, sigma = 1)

cop <- gamboostLSS:::get_copula(copula)
u1 <- gamboostLSS:::get_marginal_cdf(marginal1, y[,1], params1)
u2 <- gamboostLSS:::get_marginal_cdf(marginal2, y[,2], params2)

num_h <- (cop$pcopula(u1, u2 + eps, theta) - cop$pcopula(u1, u2 - eps, theta)) /
  (2*eps)
ana_h <- cop$h(u1, u2, theta)
stopifnot(max(abs(num_h - ana_h)) < 1e-3)

