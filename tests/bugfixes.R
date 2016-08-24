###
# Bugfixes

require("gamboostLSS")
require("gamlss")

## subset method was missing if initial mstop = 1
set.seed(1907)
x1 <- rnorm(1000)
x2 <- rnorm(1000)
x3 <- rnorm(1000)
x4 <- rnorm(1000)
x5 <- rnorm(1000)
x6 <- rnorm(1000)
mu    <- exp(1.5 + x1^2 +0.5 * x2 - 3 * sin(x3) -1 * x4)
sigma <- exp(-0.2 * x4 +0.2 * x5 +0.4 * x6)
y <- numeric(1000)
for( i in 1:1000)
    y[i] <- rnbinom(1, size = sigma[i], mu = mu[i])
dat <- data.frame(x1, x2, x3, x4, x5, x6, y)

model <- gamboostLSS(y ~ ., families = NBinomialLSS(), data = dat,
                     control = boost_control(mstop = 1))
model[10]


## selected() was broken (didn't call mboost-function)
stopifnot(all.equal(selected(model),
                    list(mu = mboost::selected(model[[1]]),
                         sigma = mboost::selected(model[[2]]))))

# If the families argument is not specified explicitly in mboostLSSone gets an
# error in cvrisk.mboostLSS() (spotted by Almond Stöcker).
# (https://github.com/boost-R/gamboostLSS/issues/9)
set.seed(1907)
x1 <- rnorm(1000)
mu    <- 2*x1
sigma <- rep(1, 1000)
y <- numeric(1000)
for( i in 1:1000)
       y[i] <- rnorm(1, mean = mu[i], sd=sigma[i])
dat <- data.frame(x1=x1, y=y)
## model with default families
model <- mboostLSS(y ~ bbs(x1), data = dat, control = boost_control(mstop = 2))
## cvrisk.mboostLSS() does not work as families was not specified in model call
cvr <- try(cvrisk(model, folds = cv(model.weights(model), B=3), trace=FALSE),
           silent = TRUE)
if (inherits(cvr, "try-error"))
    stop("cvrisk does not work if no family was (explicitly) chosen")


# Make sure that gamlss.dist::BB and gamlss.dist::BI work (spotted by F. Scheipl)
# (https://github.com/boost-R/gamboostLSS/issues/12)
set.seed(123)
n <- 100
x <- rnorm(n)
z <- rnorm(n)
data <- data.frame(y = rbinom(n, p = plogis(x + z), size = 60), x = x, z= z)
data$ymat <- with(data, cbind(success = data$y, fail = 60 - data$y))

m1 <- gamlss(ymat ~ x + z, data = data, family = BB)
m2 <- gamlss(ymat ~ x + z, data = data, family = BI)
# same with boosting
m3 <- glmboostLSS(ymat ~ x + z, data = data, families = as.families("BB"))
m4 <- glmboost(ymat ~ x + z, data = data, family = as.families("BI"))

round(data.frame(BB_gamlss = coef(m1),
                 BI_gamlss = coef(m2),
                 BB_gamboostLSS = coef(m3, off2int = TRUE, parameter = "mu"),
                 BI_gamboostLSS = coef(m4, off2int = TRUE, parameter = "mu")), 3)
