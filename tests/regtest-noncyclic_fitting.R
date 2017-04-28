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


## test that mstop = 0 is possible
compare_models <- function (m1, m2) {
  stopifnot(all.equal(coef(m1), coef(m2)))
  stopifnot(all.equal(predict(m1), predict(m2)))
  stopifnot(all.equal(fitted(m1), fitted(m2)))
  stopifnot(all.equal(selected(m1), selected(m2)))
  stopifnot(all.equal(risk(m1), risk(m2)))
  ## remove obvious differences from objects
  m1$control <- m2$control <- NULL
  m1$call <- m2$call <- NULL
  if (!all.equal(m1, m2))
    stop("Objects of offset model + 1 step and model with 1 step not identical")
  invisible(NULL)
}

# set up models
mod <- glmboostLSS(y ~ ., data = dat, method = "noncyclic", control = boost_control(mstop = 0))
mod2 <- glmboostLSS(y ~ ., data = dat, method = "noncyclic",  control = boost_control(mstop = 1))
mod3 <- glmboostLSS(y ~ ., data = dat, method = "noncyclic", control = boost_control(mstop = 1))

lapply(coef(mod), function(x) stopifnot(is.null(x)))

mstop(mod3) <- 0
mapply(compare_models, m1 = mod, m2 = mod3)

mstop(mod) <- 1
mapply(compare_models, m1 = mod, m2 = mod2)

mstop(mod3) <- 1
mapply(compare_models, m1 = mod, m2 = mod3)

## check selected
set.seed(1907)
x1 <- rnorm(500)
x2 <- rnorm(500)
x3 <- rnorm(500)
x4 <- rnorm(500)
x5 <- rnorm(500)
x6 <- rnorm(500)
mu    <- exp(1.5 +1 * x1 +0.5 * x2 -0.5 * x3 -1 * x4)
sigma <- exp(-0.4 * x3 -0.2 * x4 +0.2 * x5 +0.4 * x6)
y <- numeric(500)
for( i in 1:500)
    y[i] <- rnbinom(1, size = sigma[i], mu = mu[i])
dat <- data.frame(x1, x2, x3, x4, x5, x6, y)

model <- glmboostLSS(y ~ ., families = NBinomialLSS(), data = dat,
                     control = boost_control(mstop = 10),
                     center = TRUE, method = "cyclic")
selected(model) # ok (at least in principle)
selected(model, merge = TRUE) # ok

model <- glmboostLSS(y ~ ., families = NBinomialLSS(), data = dat,
                     control = boost_control(mstop = 10),
                     center = TRUE, method = "noncyclic")
selected(model) # ok (at least in principle)
selected(model, merge = TRUE) ## BROKEN
