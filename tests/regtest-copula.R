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

### numerical gradient check of gaussian copula 
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

### numerical gradient check of Gumbel copula 
cop <- gamboostLSS:::get_copula("gumbel")
test_points_gumbel <- list(
  # standard case - all params non-problematic
  list(u1 = 0.3, u2 = 0.6, theta = 1.4), 
  # negative correlation
  list(u1 = 0.7, u2 = 0.2, theta = 3),
  # strong dependence 
  list(u1 = 0.5, u2 = 0.5, theta = 1.5),
  # no correlation
  list(u1 = 0.1, u2 = 0.9, theta = 4)
)
eps <- 1e-5

# check analytical against numerical gradient of Gumbel copula (u1)
for (p in test_points_gumbel){
  num_grad_u1 <- (cop$logdcopula(p$u1 + eps, p$u2, p$theta) - 
                    cop$logdcopula(p$u1 - eps, p$u2, p$theta)) / (2*eps)
  
  stopifnot(abs(num_grad_u1 - cop$dlogdcopula_u1(p$u1, p$u2, p$theta)) < 1e-3)

  num_grad_theta <- (cop$logdcopula(p$u1, p$u2, p$theta + eps) -
                       cop$logdcopula(p$u1, p$u2, p$theta - eps)) / (2*eps)
  
  stopifnot(abs(num_grad_theta - cop$dlogdcopula_theta(p$u1, p$u2, p$theta)) < 1e-3)
}

### check boundary consistency of clayton copula
cop <- gamboostLSS:::get_copula("clayton")
u <- c(0.2, 0.5, 0.8)
stopifnot(max(abs(cop$pcopula(u, 1, 2) - u)) < 1e-8)
stopifnot(max(abs(cop$pcopula(1, u, 2) - u)) < 1e-8)


### check h against numeric derivative of pcopula for clayton copula
eps <- 1e-6
theta <- 2
u1 <- 0.4; u2 <- 0.6
num_h <- (cop$pcopula(u1, u2 + eps, theta) - cop$pcopula(u1, u2 - eps, theta)) / (2*eps)
ana_h <- cop$h(u1, u2, theta)
stopifnot(abs(num_h - ana_h) < 1e-5)

### check boundary consistency of gumbel copula
cop <- gamboostLSS:::get_copula("gumbel")
u <- c(0.2, 0.5, 0.8)
stopifnot(max(abs(cop$pcopula(u, 1, 2) - u)) < 1e-8)
stopifnot(max(abs(cop$pcopula(1, u, 2) - u)) < 1e-8)

### check h against numeric derivative of pcopula for gumbel copula
eps <- 1e-6
theta <- 2
u1 <- 0.4; u2 <- 0.6
num_h <- (cop$pcopula(u1, u2 + eps, theta) - cop$pcopula(u1, u2 - eps, theta)) / (2*eps)
ana_h <- cop$h(u1, u2, theta)
stopifnot(abs(num_h - ana_h) < 1e-5)

### check boundary consistency of frank copula
cop <- gamboostLSS:::get_copula("frank")
u <- c(0.2, 0.5, 0.8)
stopifnot(max(abs(cop$pcopula(u, 1, 2) - u)) < 1e-8)
stopifnot(max(abs(cop$pcopula(1, u, 2) - u)) < 1e-8)

### numerical gradient check of Frank copula
cop <- gamboostLSS:::get_copula("frank")
test_points_frank <- list(
  list(u1 = 0.3, u2 = 0.6, theta = 2),
  list(u1 = 0.7, u2 = 0.2, theta = 5),
  list(u1 = 0.5, u2 = 0.5, theta = 0.5),
  list(u1 = 0.1, u2 = 0.9, theta = 1),
  list(u1 = 0.5, u2 = 0.5, theta = -3)
)

for (p in test_points_frank){
  num_u1 <- numDeriv::grad(function(u) cop$logdcopula(u, p$u2, p$theta), p$u1)
  stopifnot(abs(num_u1 - cop$dlogdcopula_u1(p$u1, p$u2, p$theta)) < 1e-4)
  num_theta <- numDeriv::grad(function(th) cop$logdcopula(p$u1, p$u2, th), p$theta)
  stopifnot(abs(num_theta - cop$dlogdcopula_theta(p$u1, p$u2, p$theta)) < 1e-4)
}

### check h against numeric derivative of pcopula for frank copula
eps <- 1e-6
theta <- 2
u1 <- 0.4; u2 <- 0.6
num_h <- (cop$pcopula(u1, u2 + eps, theta) - cop$pcopula(u1, u2 - eps, theta)) / (2*eps)
ana_h <- cop$h(u1, u2, theta)
stopifnot(abs(num_h - ana_h) < 1e-5)

### check logdcopula against numerical mixed partial derivative of pcopula 
eps <- 1e-4
thetas <- list(gaussian = 0.5, clayton = 2, gumbel = 3, frank = 2.5)
test_u12 <- list(c(0.3,0.4), c(0.5, 0.5), c(0.7, 0.2), c(0.1, 0.9))

for (cop_name in names(thetas)){
  cop <- gamboostLSS:::get_copula(cop_name)
  theta <- thetas[[cop_name]]
  for (p in test_u12){
    u1 <- p[1]; u2 <- p[2]
    num_c <- (cop$pcopula(u1+eps, u2+eps, theta) - 
                cop$pcopula(u1+eps, u2-eps, theta) -
                cop$pcopula(u1-eps, u2+eps, theta) +
                cop$pcopula(u1-eps, u2-eps, theta)) / (4*eps^2)
    ana_logc <- cop$logdcopula(u1, u2, theta)
    stopifnot(abs(num_c - exp(ana_logc)) < 1e-4)
  }
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

### check h against numeric derivative of pcopula for gaussian copula
set.seed(42)
y <- cbind(rnorm(10), rnorm(10))
eps <- 1e-3
theta <- 0.5
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


### check binary_continuous_loglik at theta = 0
y_bc <- cbind(c(0, 1, 0, 1), rnorm(4))
mu1_val <- 0.3
params_bin <- list(mu = mu1_val)
params_norm <- list(mu = 0, sigma = 1)

loglik_bc <- gamboostLSS:::binary_continuous_loglik(
  y_bc, BernoulliLSS(), GaussianLSS(), params_bin, params_norm, "gaussian", 0)

expected_bc <- log(ifelse(y_bc[,1] == 1, mu1_val, 1 - mu1_val)) + 
  dnorm(y_bc[,2], log = TRUE)

stopifnot(max(abs(loglik_bc - expected_bc)) < 1e-6)


### reference value test for binary_continuous_loglik
y_ref <- cbind(0, 1)
ref_val <- gamboostLSS:::binary_continuous_loglik(
  y_ref, BernoulliLSS(), GaussianLSS(), list(mu = 0.4), list(mu = 0, sigma =1),
  "gaussian", 0.5)

stopifnot(abs(ref_val - (-2.36596)) < 1e-2)

  
### test correct naming of sub-family parameters
stopifnot(identical(names(gamboostLSS:::CopulaFamilies(GaussianLSS(),
                                                       GaussianLSS())),
                    c("mu1", "sigma1", "mu2", "sigma2", "theta")))

stopifnot(identical(names(gamboostLSS:::CopulaFamilies(NBinomialLSS(),
                                                       NBinomialLSS())),
                    c("mu1", "sigma1", "mu2", "sigma2", "theta")))

stopifnot(identical(names(gamboostLSS:::CopulaFamilies(
  gamboostLSS:::BernoulliLSS(), 
  gamboostLSS:::BernoulliLSS())),
  c("mu1", "mu2", "theta")))

stopifnot(identical(names(gamboostLSS:::CopulaFamilies(
  gamboostLSS:::BernoulliLSS(), GaussianLSS())),
  c("mu1", "mu2", "sigma2", "theta")))

stopifnot(identical(names(gamboostLSS:::CopulaFamilies(
  GaussianLSS(), gamboostLSS:::BernoulliLSS())),
  c("mu1", "sigma1", "mu2", "theta")))


### check CopulaFamilies copula response/name
for (cop_name in c("gaussian", "clayton", "gumbel", "frank")){
  copula_family <- gamboostLSS:::CopulaFamilies(GaussianLSS(), GaussianLSS(),
                                                copula = cop_name)
  cop <- gamboostLSS:::get_copula(cop_name)
  stopifnot(identical(copula_family$theta@name, cop$name))
  
  # compare responses
  f_vals <- c(-2, -0.5, 0, 0.5, 2)
  stopifnot(max(abs(copula_family$theta@response(f_vals) - 
                      cop$response(f_vals))) < 1e-5)
}

### check "unknown-copula" case for get_copula
result <- tryCatch(gamboostLSS:::get_copula("unknown"),
                   error = function(e) "error")
stopifnot(identical(result, "error"))


### check ngradient_theta against numDeriv for gaussian copula
library(numDeriv)
set.seed(42)
y <- cbind(rnorm(10), rnorm(10))
cf <- gamboostLSS:::CopulaFamilies(GaussianLSS(), GaussianLSS(), 
                                   copula = "gaussian")

w <- rep(1, nrow(y))
invisible(cf$mu1@offset(y, w))
invisible(cf$sigma1@offset(y, w))
invisible(cf$mu2@offset(y, w))
invisible(cf$sigma2@offset(y, w))

f0 <- 0.3
num_grad <- numDeriv::grad(function(f) -cf$theta@risk(y, f), f0)
ana_grad <- sum(cf$theta@ngradient(y, f0))
stopifnot(abs(num_grad - ana_grad) < 1e-5)

### check ngradient_theta against numDeriv for clayton copula
set.seed(42)
y <- cbind(rnorm(10), rnorm(10))
cf <- gamboostLSS:::CopulaFamilies(GaussianLSS(), GaussianLSS(), 
                                   copula = "clayton")

w <- rep(1, nrow(y))
invisible(cf$mu1@offset(y, w))
invisible(cf$sigma1@offset(y, w))
invisible(cf$mu2@offset(y, w))
invisible(cf$sigma2@offset(y, w))

f0 <- 0.5
num_grad <- numDeriv::grad(function(f) - cf$theta@risk(y, f), f0)
ana_grad <- sum(cf$theta@ngradient(y, f0))
stopifnot(abs(num_grad - ana_grad) < 1e-4)

### check ngradient_theta against numDeriv for gumbel copula
set.seed(42)
y <- cbind(rnorm(10), rnorm(10))
cf <- gamboostLSS:::CopulaFamilies(GaussianLSS(), GaussianLSS(), 
                                   copula = "gumbel")

w <- rep(1, nrow(y))
invisible(cf$mu1@offset(y, w))
invisible(cf$sigma1@offset(y, w))
invisible(cf$mu2@offset(y, w))
invisible(cf$sigma2@offset(y, w))

f0 <- 0.5
num_grad <- numDeriv::grad(function(f) - cf$theta@risk(y, f), f0)
ana_grad <- sum(cf$theta@ngradient(y, f0))
stopifnot(abs(num_grad - ana_grad) < 1e-4)

### check get_marginal_dcdf analytical formula for gaussian
eps <- 1e-6
y <- c(0.5, -1.2, 2.0)
params <- list(mu = 0.3, sigma = 1.5)

num_dmu <- (gamboostLSS:::get_marginal_cdf(GaussianLSS(), y, 
                                            list(mu = params$mu + eps,
                                                 sigma = params$sigma)) -
              gamboostLSS:::get_marginal_cdf(GaussianLSS(), y, 
                                              list(mu = params$mu - eps,
                                                   sigma = params$sigma))) / 
  (2*eps)
ana_dmu <- gamboostLSS:::get_marginal_dcdf(GaussianLSS(), y, params, "mu")
stopifnot(max(abs(num_dmu - ana_dmu)) < 1e-5)

num_dsigma <- (gamboostLSS:::get_marginal_cdf(GaussianLSS(), y, 
                                              list(mu = params$mu,
                                                   sigma = params$sigma + eps)) -
                 gamboostLSS:::get_marginal_cdf(GaussianLSS(), y, 
                                                list(mu = params$mu,
                                                     sigma = params$sigma- eps))) / 
  (2*eps)
ana_dsigma <- gamboostLSS:::get_marginal_dcdf(GaussianLSS(), y, params, "sigma")
stopifnot(max(abs(num_dsigma - ana_dsigma)) < 1e-5)

### check get_marginal_dcdf fallback for negative binomial
y_nb <- c(1, 3, 5)
params_nb <- list(mu = 2, sigma = 1)
d_nb <- gamboostLSS:::get_marginal_dcdf(NBinomialLSS(), y_nb, params_nb, "mu")
stopifnot(all(is.finite(d_nb)))

### check ngradient_p (mu1, sigma1) against numDeriv for gaussian copula 
set.seed(42)
y <- cbind(rnorm(10, 0.5, 1.5), rnorm(10))
cf <- gamboostLSS:::CopulaFamilies(GaussianLSS(), GaussianLSS(),
                                   copula = "gaussian", theta = 0)

w <- rep(1, nrow(y))
invisible(cf$mu1@offset(y,w))
invisible(cf$sigma1@offset(y,w))
invisible(cf$mu2@offset(y,w))
invisible(cf$sigma2@offset(y,w))

f0 <- 0.3
num_grad_mu1 <- numDeriv::grad(function(f) -cf$mu1@risk(y, f, w=w), f0)
ana_grad_mu1 <- sum(cf$mu1@ngradient(y, f0, w))
stopifnot(abs(num_grad_mu1 - ana_grad_mu1) < 1e-4)

f0_sigma <- log(1.5)
num_grad_sigma1 <- numDeriv::grad(function(f) - cf$sigma1@risk(y, f, w=w),
                                  f0_sigma)
ana_grad_sigma1 <- sum(cf$sigma1@ngradient(y, f0_sigma, w))
stopifnot(abs(num_grad_sigma1 - ana_grad_sigma1) < 1e-4)

### check ngradient_p for mu2 and sigma2 against numDeriv
num_grad_mu2 <- numDeriv::grad(function(f) -cf$mu2@risk(y, f, w = w), 0.3)
ana_grad_mu2 <- sum(cf$mu2@ngradient(y, 0.3, w))
stopifnot(abs(num_grad_mu2 - ana_grad_mu2) < 1e-4)

num_grad_sigma2 <- numDeriv::grad(function(f) -cf$sigma2@risk(y, f, w = w), 
                                  log(1.5))
ana_grad_sigma2 <- sum(cf$sigma2@ngradient(y, log(1.5), 2))
stopifnot(abs(num_grad_sigma2 - ana_grad_sigma2) < 1e-4)

### verify all ngradients return finite values after initialization
set.seed(42)
y <- cbind(rnorm(20), rnorm(20))
cf <- gamboostLSS:::CopulaFamilies(GaussianLSS(), GaussianLSS(), 
                                   copula = "gaussian", theta = 0)
w <- rep(1, nrow(y))

invisible(cf$mu1@offset(y, w))
invisible(cf$sigma1@offset(y, w))
invisible(cf$mu2@offset(y, w))
invisible(cf$sigma2@offset(y, w))
invisible(cf$theta@offset(y, w))

stopifnot(all(is.finite(cf$mu1@ngradient(y, 0.5, w))))
stopifnot(all(is.finite(cf$sigma1@ngradient(y, log(1.5), w))))
stopifnot(all(is.finite(cf$mu2@ngradient(y, 0.5, w))))
stopifnot(all(is.finite(cf$sigma2@ngradient(y, log(1.5), w))))
stopifnot(all(is.finite(cf$theta@ngradient(y, 0.3, w))))

### check ngradient_theta against numDeriv for Frank copula
set.seed(42)
y <- cbind(rnorm(10), rnorm(10))
cf <- gamboostLSS:::CopulaFamilies(GaussianLSS(), GaussianLSS(),
                                   copula = "frank")

w <- rep(1, nrow(y))
invisible(cf$mu1@offset(y, w))
invisible(cf$sigma1@offset(y, w))
invisible(cf$mu2@offset(y, w))
invisible(cf$sigma2@offset(y, w))

f0 <- 2
num_grad <- numDeriv::grad(function(f) - cf$theta@risk(y, f), f0)
ana_grad <- sum(cf$theta@ngradient(y, f0))
stopifnot(abs(num_grad - ana_grad) < 1e-4)

### check theta offsets
set.seed(42)
y <- cbind(rnorm(40), rnorm(40))
y[, 2] <- y[, 1] + rnorm(40)
w <- rep(1, 40)
tau_hat <- cor(y[, 1], y[, 2], method = "kendall")

cf <- gamboostLSS:::CopulaFamilies(GaussianLSS(), GaussianLSS(),
                                   copula = "gaussian")
stopifnot(abs(cf$theta@offset(y, w) - atanh(sin(pi * tau_hat / 2))) < 1e-10)

cf <- gamboostLSS:::CopulaFamilies(GaussianLSS(), GaussianLSS(),
                                   copula = "clayton")
stopifnot(abs(cf$theta@offset(y, w) - log(2 * tau_hat / (1 - tau_hat))) < 1e-10)

cf <- gamboostLSS:::CopulaFamilies(GaussianLSS(), GaussianLSS(),
                                   copula = "gumbel")
stopifnot(abs(cf$theta@offset(y, w) - log(tau_hat / (1 - tau_hat))) < 1e-10)

cf <- gamboostLSS:::CopulaFamilies(GaussianLSS(), GaussianLSS(),
                                   copula = "frank")
stopifnot(abs(cop$tau(cf$theta@offset(y, w)) - tau_hat) < 1e-6)

### check fallback for negative dependence in Clayton/Gumbel 
y_neg <- y; y_neg[, 2] <- -y[, 2]
cf <- gamboostLSS:::CopulaFamilies(GaussianLSS(), GaussianLSS(),
                                  copula = "clayton")
stopifnot(cf$theta@offset(y_neg, w) == 0)

cf <- gamboostLSS:::CopulaFamilies(GaussianLSS(), GaussianLSS(),
                                  copula = "gumbel")
stopifnot(cf$theta@offset(y_neg, w) == 0)


### end-to-end test 
set.seed(42)
n <- 100
x <- rnorm(n)
rho <- 0.5
eps <- matrix(rnorm(2*n), n, 2) %*% chol(matrix(c(1, rho, rho, 1), 2, 2))
y <- cbind(x + eps[, 1], -x + eps[,2])
df <- data.frame(x = x)

for (cop_name in c("gaussian", "clayton", "gumbel", "frank")){
  cf <- gamboostLSS:::CopulaFamilies(GaussianLSS(), GaussianLSS(), 
                                     copula = cop_name)
  fit <- gamboostLSS(y ~ x, families = cf, data = df, 
                     control = boost_control(mstop = 10, nu = 0.1))
  stopifnot(all(sapply(fit, function(m) all(is.finite(fitted(m))))))
}

### check stabilization methods for theta ngradient  
set.seed(42)
y <- cbind(rnorm(20), rnorm(20))
w <- rep(1, 20)

make_cf <- function(stabilization){
  cf <- gamboostLSS:::CopulaFamilies(GaussianLSS(), GaussianLSS(), 
                                     copula = "gaussian", 
                                     stabilization = stabilization)
  for (p in c("mu1", "sigma1", "mu2", "sigma2", "theta"))
    invisible(cf[[p]]@offset(y, w))
  cf
}

cf_none <- make_cf("none")
cf_mad <- make_cf("MAD")
cf_l2 <- make_cf("L2")

f0 <- 0.3
raw <- cf_none$theta@ngradient(y, f0, w)

div_mad <- weighted.median(
  abs(raw - weighted.median(raw, w = w)), w = w)
stopifnot(all.equal(cf_mad$theta@ngradient(y, f0, w), raw/div_mad))

div_l2 <- sqrt(weighted.mean(raw^2, w = w))
stopifnot(all.equal(cf_l2$theta@ngradient(y, f0, w), raw / div_l2))

### check stabilization for marginal ngradient
raw_mu <- cf_none$mu1@ngradient(y, f0, w)
div_mu <- weighted.median(abs(raw_mu - weighted.median(raw_mu, w = w)), w = w)
stopifnot(all.equal(cf_mad$mu1@ngradient(y, f0, w), raw_mu / div_mu))


### check Kendall's tau conversion
cop <- gamboostLSS:::get_copula("gaussian")
stopifnot(abs(cop$tau(0.5) - 2/pi * asin(0.5)) < 1e-10)
for (th in c(-0.8, -0.3, 0.2, 0.9))
  stopifnot(abs(cop$tau_inv(cop$tau(th)) - th) < 1e-8)

cop <- gamboostLSS:::get_copula("clayton")
stopifnot(abs(cop$tau(2) - 0.5) < 1e-10)
for (th in c(0.5, 1, 2, 8))
  stopifnot(abs(cop$tau_inv(cop$tau(th)) - th) < 1e-8)

cop <- gamboostLSS:::get_copula("gumbel")
stopifnot(abs(cop$tau(2) - 0.5) < 1e-10)
for (th in c(1.2, 2, 5, 10))
  stopifnot(abs(cop$tau_inv(cop$tau(th)) - th) < 1e-8)

cop <- gamboostLSS:::get_copula("frank")
stopifnot(abs(cop$tau(3) + cop$tau(-3)) < 1e-8)
stopifnot(cop$tau(8) > cop$tau(2))
for (th in c(-8, -2, 1, 5, 20))
  stopifnot(abs(cop$tau_inv(cop$tau(th)) - th) < 1e-6)

### check response_inv
for (cn in c("gaussian", "clayton", "gumbel", "frank")){
  cop <- gamboostLSS:::get_copula(cn)
  for (f in c(-1, 0.3, 2))
    stopifnot(abs(cop$response_inv(cop$response(f)) - f) < 1e-10)
}




