### 
# Test copula kernel functions: logdcopula, dlogdcopula_u1, dlogdcopula_theta

require("gamboostLSS")
library(copula)

u1 <- c(0.2, 0.5, 0.8)
u2 <- c(0.3, 0.5, 0.7)

### check gaussian copula density against copula-package 
# check for correlation parameter rho = 0.5
ref <- dCopula(cbind(u1, u2), normalCopula(0.5), log = TRUE)
own <- gamboostLSS:::get_copula("gaussian")$logdcopula(u1, u2, 0.5)
stopifnot(max(abs(ref - own)) < 1e-8)

### check clayton copula density against copula-package
# check for correlation parameter theta = 2
ref <- dCopula(cbind(u1, u2), claytonCopula(2), log = TRUE)
own <- gamboostLSS:::get_copula("clayton")$logdcopula(u1, u2, 2)
stopifnot(max(abs(ref - own)) < 1e-8)

### check gumbel copula density against copula-package
# check for correlation parameter theta = 2
ref <- dCopula(cbind(u1, u2), gumbelCopula(2), log = TRUE)
own <- gamboostLSS:::get_copula("gumbel")$logdcopula(u1, u2, 2)
stopifnot(max(abs(ref - own)) < 1e-8)

### check frank copula density against copula-package
# check for correlation parameter theta = 2
ref <- dCopula(cbind(u1, u2), frankCopula(2), log = TRUE)
own <- gamboostLSS:::get_copula("frank")$logdcopula(u1, u2, 2)
stopifnot(max(abs(ref - own)) < 1e-8)

# ------------------------------------------------------------------------------

eps <- 1e-6
### check numerical gradient stability for gaussian copula
cop <- gamboostLSS:::get_copula("gaussian")

num_grad <- (cop$logdcopula(0.3 + eps, 0.6, 0.4) - 
               cop$logdcopula(0.3 - eps, 0.6, 0.4)) / (2*eps)
stopifnot(abs(num_grad - cop$dlogdcopula_u1(0.3, 0.6, 0.4)) < 1e-5)

### check numerical gradient stability for clayton copula
cop <- gamboostLSS:::get_copula("clayton")

num_grad <- (cop$logdcopula(0.3 + eps, 0.6, 0.4) -
               cop$logdcopula(0.3 - eps, 0.6, 0.4)) / (2*eps)
stopifnot(abs(num_grad - cop$dlogdcopula_u1(0.3, 0.6, 0.4)) < 1e-5)

### check numerical gradient stability for gumbel copula
cop <- gamboostLSS:::get_copula("gumbel")

num_grad <- (cop$logdcopula(0.3 + eps, 0.6, 2) -
               cop$logdcopula(0.3 - eps, 0.6, 2)) / (2*eps)
stopifnot(abs(num_grad - cop$dlogdcopula_u1(0.3, 0.6, 2)) < 1e-5)

### check numerical gradient stability for frank copula
cop <- gamboostLSS:::get_copula("frank")

num_grad <- (cop$logdcopula(0.3 + eps, 0.6, 0.4) -
               cop$logdcopula(0.3 - eps, 0.6, 0.4)) / (2*eps)
stopifnot(abs(num_grad - cop$dlogdcopula_u1(0.3, 0.6, 0.4)) < 1e-5)

# ------------------------------------------------------------------------------

### check numerical stability at singularities
for (nm in c("gaussian", "clayton", "gumbel", "frank")){
  cop <- gamboostLSS:::get_copula(nm)
  val <- cop$logdcopula(1e-6, 0.5, ifelse(nm == "gaussian", 0.5, 2))
  stopifnot(is.finite(val))
}

# ------------------------------------------------------------------------------

### check input stability
tryCatch(gamboostLSS:::get_copula("xy"),
         error = function(e) cat("unknown copula throws error\n"))

# ------------------------------------------------------------------------------

### check get_marginal_cdf
fam <- GaussianLSS()
u <- gamboostLSS:::get_marginal_cdf(fam, y = c(-1, 0, 1), params = list(mu = 0, sigma = 1))

# domain restrictions for cdf
stopifnot(all(is.finite(u)), all(u > 0), all(u < 1))

# deviance from predefined cdf must be marginal (no pun intended)
stopifnot(max(abs(u - pnorm(c(-1, 0, 1)))) < 1e-10)

# ------------------------------------------------------------------------------

### check CopulaFamilies structure
fam <- CopulaFamilies(GaussianLSS(), GaussianLSS(), copula = "gaussian",
                      theta = 0.5)
# class
stopifnot(inherits(fam, "families"))

# slots
stopifnot(all(c("mu1", "sigma1", "mu2", "sigma2", "theta_cop") %in% 
                names(fam)))

# slots are valid objects
for (nm in names(fam))
  stopifnot(inherits(fam[[nm]], "boost_family"))

# check_y takes only matrix
tryCatch(fam$mu1@check_y(c(1, 2, 3)),
         error = function(e) cat("check_y non-matrix: OK\n"))

# check_y accepts matrix
stopifnot(is.matrix(fam$mu1@check_y(cbind(c(1, 2), c(3, 4)))))


