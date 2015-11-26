require("gamboostLSS")

###negbin dist, linear###


#there was a bug that if mstop is length(families) + 1
set.seed(2611)
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

model <- glmboostLSS(y ~ ., families = NBinomialLSS(), data = dat,
                     control = boost_control(mstop = 3), cycling = FALSE)


set.seed(2611)
#Regular fit and fit with risk = "oob" should be identical
w <- sample(c(0,1), 1000, replace = TRUE)

model <- glmboostLSS(y ~ ., families = NBinomialLSS(), data = dat,
                     control = boost_control(mstop = 100), cycling = TRUE,
                     weights = w)

model_oob <- glmboostLSS(y ~ ., families = NBinomialLSS(), data = dat,
                         control = boost_control(mstop = 100, risk = "oob"), 
                         cycling = TRUE,
                         weights = w)

stopifnot(all(unlist(model$mu$coef()) ==  unlist(model_oob$mu$coef())) &
          all(unlist(model$sigma$coef()) ==  unlist(model_oob$sigma$coef())))
