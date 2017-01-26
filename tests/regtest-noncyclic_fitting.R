require("gamboostLSS")

###negbin dist, linear###

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
for (i in 1:1000)
  y[i] <- rnbinom(1, size = sigma[i], mu = mu[i])
dat <- data.frame(x1, x2, x3, x4, x5, x6, y)

#fit models at number of params + 1

#glmboost
model <- glmboostLSS(y ~ ., families = NBinomialLSS(), data = dat,
                     control = boost_control(mstop = 3), method = "noncyclic")

#linear baselearner with bols
model <- gamboostLSS(y ~ ., families = NBinomialLSS(), data = dat,
                     control = boost_control(mstop = 3), method = "noncyclic",
                     baselearner = "bols")

#nonlinear bbs baselearner

model <- gamboostLSS(y ~ ., families = NBinomialLSS(), data = dat,
                     control = boost_control(mstop = 3), method = "noncyclic",
                     baselearner = "bbs")

#reducing model and increasing it afterwards should yield the same fit

model <- glmboostLSS(y ~ ., families = NBinomialLSS(), data = dat,
                     control = boost_control(mstop = 50), method = "noncyclic")

m_co <- coef(model)

mstop(model) <- 5
mstop(model) <- 50

stopifnot(all.equal(m_co, coef(model)))


model <- gamboostLSS(y ~ ., families = NBinomialLSS(), data = dat,
                     control = boost_control(mstop = 50), method = "noncyclic",
                     baselearner = "bols")

m_co <- coef(model)

mstop(model) <- 5
mstop(model) <- 50

stopifnot(all.equal(m_co, coef(model)))


model <- gamboostLSS(y ~ ., families = NBinomialLSS(), data = dat,
                     control = boost_control(mstop = 50), method = "noncyclic",
                     baselearner = "bbs")

m_co <- coef(model)

mstop(model) <- 5
mstop(model) <- 50

stopifnot(all.equal(m_co, coef(model)))


model <- gamboostLSS(y ~ ., families = NBinomialLSS(), data = dat,
                     control = boost_control(mstop = 50), method = "noncyclic",
                     baselearner = "bbs")

m_co <- coef(model)

mstop(model) <- 5
mstop(model) <- 50

stopifnot(all.equal(m_co, coef(model)))

## check cvrisk for noncyclic models
model <- glmboostLSS(y ~ ., families = NBinomialLSS(), data = dat,
                     control = boost_control(mstop = 3), method = "noncyclic")
cvr1 <- cvrisk(model, grid = 1:50, cv(model.weights(model), B = 5))
cvr1
plot(cvr1)

risk(model, merge = TRUE)
risk(model, merge = FALSE)