###
# check glmboostLSS()

require("gamboostLSS")
require("gamlss")

set.seed(1907)
n <- 5000
x1  <- runif(n)
x2 <- runif(n)
mu <- 2 -1*x1 - 3*x2
sigma <- exp(-1*x1 + 3*x2)
df <- exp(1 + 3*x1 + 1*x2)
y <- rTF(n = n, mu = mu, sigma = sigma, nu = df)

### check subset method
model <- glmboostLSS(y ~ x1 + x2, families = StudentTLSS(),
                     control = boost_control(mstop = 10),
                     center = TRUE)
model2 <- glmboostLSS(y ~ x1 + x2, families = StudentTLSS(),
                          control = boost_control(mstop = 20),
                          center = TRUE)
model[20]
stopifnot(all.equal(coef(model),coef(model2)))
stopifnot(length(coef(model2, aggregate = "none")[[1]][[1]]) ==
          length(coef(model, aggregate = "none")[[1]][[1]]))
stopifnot(length(coef(model2, aggregate = "none")[[1]][[1]]) == 20)

model <- glmboostLSS(y ~ x1 + x2, families = StudentTLSS(),
                     control = boost_control(mstop = 10),
                     center = TRUE)
model2[10]
stopifnot(all.equal(coef(model),coef(model2)))
stopifnot(length(coef(model2, aggregate = "none")[[1]][[1]]) ==
          length(coef(model, aggregate = "none")[[1]][[1]]))
stopifnot(length(coef(model2, aggregate = "none")[[1]][[1]]) == 10)

### check trace
model <- glmboostLSS(y ~ x1 + x2, families = StudentTLSS(),
                     control = boost_control(mstop = 10, trace =TRUE),
                     center = TRUE)
model[20]

### check formula-interface with lists
set.seed(1907)
n <- 5000
x1  <- runif(n)
x2 <- runif(n)
mu <- 2 - 3*x2
sigma <- exp(-1*x1 + 3*x2)
df <- exp(1 + 3*x1)
y <- rTF(n = n, mu = mu, sigma = sigma, nu = df)
model <- glmboostLSS(list(mu = y ~ x2,
                          sigma = y ~ x1 + x2,
                          df = y ~ x1),
                     families = StudentTLSS(),
                     control = boost_control(mstop = 10, trace =TRUE),
                     center = TRUE)

stopifnot(all.equal(lapply(coef(model, which = ""), function(x) names(x)[-1]),
                    list(mu = "x2", sigma = c("x1", "x2"), df = "x1")))

model <- glmboostLSS(list(mu = y ~ x2,
                          df = y ~ x1,
                          sigma = y ~ x1 + x2),
                     families = StudentTLSS(),
                     control = boost_control(mstop = 10, trace =TRUE),
                     center = TRUE)

stopifnot(all.equal(lapply(coef(model, which = ""), function(x) names(x)[-1]),
                    list(mu = "x2", sigma = c("x1", "x2"), df = "x1")))


### check weights interface
set.seed(1907)
n <- 2500
x1  <- runif(n)
x2 <- runif(n)
mu <- 2 - 3*x2
sigma <- exp(-1*x1 + 3*x2)
df <- exp(1 + 3*x1)
y <- rTF(n = n, mu = mu, sigma = sigma, nu = df)
dat <- data.frame(x1, x2, y)
dat2 <- rbind(dat, dat) # data frame with duplicate entries

model <- glmboostLSS(list(mu = y ~ x2,
                          sigma = y ~ x1 + x2,
                          df = y ~ x1),
                     data = dat, weights = rep(2, nrow(dat)),
                     families = StudentTLSS(),
                     control = boost_control(mstop = 10, trace =TRUE),
                     center = TRUE)

model2 <- glmboostLSS(list(mu = y ~ x2,
                           sigma = y ~ x1 + x2,
                           df = y ~ x1),
                      data = dat2, families = StudentTLSS(),
                      control = boost_control(mstop = 10, trace =TRUE),
                      center = TRUE)

stopifnot(all.equal(coef(model),coef(model2)))
