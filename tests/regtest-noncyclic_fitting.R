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
                     control = boost_control(mstop = 3), method = "outer")

model <- glmboostLSS(y ~ ., families = NBinomialLSS(), data = dat,
                     control = boost_control(mstop = 3), method = "inner")

#linear baselearner with bols
model <- gamboostLSS(y ~ ., families = NBinomialLSS(), data = dat,
                     control = boost_control(mstop = 3), method = "outer",
                     baselearner = "bols")

model <- gamboostLSS(y ~ ., families = NBinomialLSS(), data = dat,
                     control = boost_control(mstop = 3), method = "inner",
                     baselearner = "bols")

#nonlinear bbs baselearner

model <- gamboostLSS(y ~ ., families = NBinomialLSS(), data = dat,
                     control = boost_control(mstop = 3), method = "outer",
                     baselearner = "bbs")

model <- gamboostLSS(y ~ ., families = NBinomialLSS(), data = dat,
                     control = boost_control(mstop = 3), method = "inner",
                     baselearner = "bbs")



#reducing model and increasing it afterwards should yield the same fit

model <- glmboostLSS(y ~ ., families = NBinomialLSS(), data = dat,
                     control = boost_control(mstop = 50), method = "outer")

m_co <- coef(model)

mstop(model) <- 5
mstop(model) <- 50

stopifnot(all.equal(m_co, coef(model)))


model <- gamboostLSS(y ~ ., families = NBinomialLSS(), data = dat,
                     control = boost_control(mstop = 50), method = "outer",
                     baselearner = "bols")

m_co <- coef(model)

mstop(model) <- 5
mstop(model) <- 50

stopifnot(all.equal(m_co, coef(model)))


model <- gamboostLSS(y ~ ., families = NBinomialLSS(), data = dat,
                     control = boost_control(mstop = 50), method = "outer",
                     baselearner = "bbs")

m_co <- coef(model)

mstop(model) <- 5
mstop(model) <- 50

stopifnot(all.equal(m_co, coef(model)))


##inner###

model <- glmboostLSS(y ~ ., families = NBinomialLSS(), data = dat,
                     control = boost_control(mstop = 50), method = "inner")

m_co <- coef(model)

mstop(model) <- 5
mstop(model) <- 50

stopifnot(all.equal(m_co, coef(model)))

model <- gamboostLSS(y ~ ., families = NBinomialLSS(), data = dat,
                     control = boost_control(mstop = 50), method = "inner",
                     baselearner = "bols")

m_co <- coef(model)

mstop(model) <- 5
mstop(model) <- 50

stopifnot(all.equal(m_co, coef(model)))


model <- gamboostLSS(y ~ ., families = NBinomialLSS(), data = dat,
                     control = boost_control(mstop = 50), method = "inner",
                     baselearner = "bbs")

m_co <- coef(model)

mstop(model) <- 5
mstop(model) <- 50

stopifnot(all.equal(m_co, coef(model)))

## check cvrisk for noncyclic models
model <- glmboostLSS(y ~ ., families = NBinomialLSS(), data = dat,
                     control = boost_control(mstop = 3), method = "outer")
cvr1 <- cvrisk(model, grid = 1:50, cv(model.weights(model), B = 5))
cvr1
plot(cvr1)

risk(model, merge = TRUE)
risk(model, merge = FALSE)

model <- glmboostLSS(y ~ ., families = NBinomialLSS(), data = dat,
                     control = boost_control(mstop = 3), method = "inner")
cvr2 <- cvrisk(model, grid = 1:50, cv(model.weights(model), B = 5))
cvr2
plot(cvr2)

risk(model, merge = TRUE)
risk(model, merge = FALSE)

model <- glmboostLSS(y ~ ., families = NBinomialLSS(), data = dat,
                     control = boost_control(mstop = 3), method = "cycling")
cvr3 <- cvrisk(model, grid = make.grid(c(mu = 25, sigma = 25)), cv(model.weights(model), B = 5))
cvr3
plot(cvr3)
plot(cvr3, type = "lines")

risk(model, merge = TRUE)
risk(model, merge = FALSE)

par(mfrow = c(2, 2))
plot(cvr1)
plot(cvr2)
plot(cvr3, type = "heatmap")
plot(cvr3, type = "lines")


