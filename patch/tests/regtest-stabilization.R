###
# test functionality of stabilization

require("gamboostLSS")

## simulate Gaussian data
set.seed(0804)
x1 <- runif(1000)
x2 <- runif(1000)
x3 <- runif(1000)
x4 <- runif(1000)
mu    <- 1.5 +1 * x1 +4 * x2
sigma <- exp(1 - 0.2 * x3 -0.4 * x4)
y <- rnorm(mean = mu, sd = sigma, n = length(mu))
dat <- data.frame(x1, x2, x3, x4, y)

## do not use stabilization
m1 <- glmboostLSS(y ~ x1 + x2 + x3 + x4,
                  families=GaussianLSS(),
                  data=dat)
coef(m1)

## use stabilization via options (for backwards compatibility)
options(gamboostLSS_stab_ngrad = TRUE)
m2 <- glmboostLSS(y ~ x1 + x2 + x3 + x4,
                  families=GaussianLSS(),
                  data=dat)
coef(m2)
options(gamboostLSS_stab_ngrad = FALSE)

## now use novel interface via families:
m3 <- glmboostLSS(y ~ x1 + x2 + x3 + x4,
                  families = GaussianLSS(stabilization = "MAD"),
                  data=dat)
stopifnot(all.equal(coef(m3), coef(m2)))

## check if everything is handled correctly
res <- try(GaussianLSS(stabilization = "test"),silent = TRUE)
res

## continue these checks for (all?) other families
