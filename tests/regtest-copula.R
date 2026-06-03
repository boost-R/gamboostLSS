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

### check analytical against numerical gradient of gaussian-copula (u1)
for (p in test_points){
  num_grad_u1 <- (cop$logdcopula(p$u1 + eps, p$u2, p$theta) - 
                    cop$logdcopula(p$u1 - eps, p$u2, p$theta)) / (2*eps)
  ana_grad_u1 <- cop$dlogdcopula_u1(p$u1, p$u2, p$theta)
  
  stopifnot(abs(num_grad_u1 - ana_grad_u1) < 1e-5)
}

### check analytical against numerical gradient of gaussian-copula (theta)
for (p in test_points){
  num_grad_theta <- (cop$logdcopula(p$u1, p$u2, p$theta + eps) -
                       cop$logdcopula(p$u1, p$u2, p$theta - eps)) / (2*eps)
  ana_grad_theta <- cop$dlogdcopula_theta(p$u1, p$u2, p$theta)
  
  stopifnot(abs(num_grad_theta - ana_grad_theta) < 1e-5)
}

