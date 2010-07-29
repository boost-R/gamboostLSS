\name{mboostLSS}
\alias{mboostLSS}
\alias{blackboostLSS}
\alias{glmboostLSS}
\alias{gamboostLSS}
\alias{mboostLSS_fit}

\title{
  Fitting mboostLSS models
}
\description{
  Functions for fitting LSS (location, shape and scale) models using
  boosting techniques.
}
\usage{
mboostLSS(formula, data = list(), families = list(),
          control = boost_control(), weights = NULL, ...)
glmboostLSS(formula, data = list(), families = list(),
            control = boost_control(), weights = NULL, ...)
gamboostLSS(formula, data = list(), families = list(),
            control = boost_control(), weights = NULL, ...)
blackboostLSS(formula, data = list(), families = list(),
              control = boost_control(), weights = NULL, ...)

## fit function:
mboostLSS_fit(formula, data = list(), families = list(),
              control = boost_control(), weights = NULL, fun = mboost, ...)
}

\arguments{
  \item{formula}{ a symbolic description of the model to be fit. See
    \code{\link{mboost}} for details. If \code{formula} is a single formula,
    the same formula is used for all model components. \code{formula}
    can also be a (named) list, where each list element corresponds
    to one component of the LSS model. The names must be the same as in
    the family (see example for details). }
  \item{data}{ a data frame containing the variables in the model.}
  \item{families}{ a list of mboost-\code{\link{Family}}s constructed
    via \code{\link{makeFamilies}}. }
  \item{control}{ a list of parameters controlling the algorithm. For
    more details see \code{\link{boost_control}}. }
  \item{weights}{ a numeric vector of weights (optional). }
  \item{fun}{ fit function. Either \code{\link{mboost}},
    \code{\link{glmboost}}, \code{\link{gamboost}} or
    \code{\link{blackboost}}. Specified directly via the according LSS
    function. E.g. \code{gamboostLSS()} calls
    \code{mboostLSS_fit(..., fun = gamboost)}. }
  \item{\dots}{ Further arguments to be passed to \code{mboostLSS_fit}.
    See there for possible arguments and details. In
    \code{mboostLSS_fit}, \code{\dots} represent further arguments to be
    passed to \code{\link{mboost}} and \code{\link{mboost_fit}}. So
    \code{\dots} can be all arguments of \code{mboostLSS_fit} and
    \code{mboost_fit}. }
}
\details{
  \code{glmboostLSS} uses \code{\link{glmboost}} to fit the components
  of an LSS (location, shape and scale) model, i.e., a component-wise
  linear boosting model is fitted for each component. Fitting is achieved
  by iterating the updates of the single components and using the fitted
  values in the other components.

  \code{gamboostLSS} uses \code{\link{gamboost}} to fit the components
  of an LSS (location, shape and scale) model, i.e., a component-wise
  additive boosting model is fitted for each component.

  \code{mboostLSS} uses \code{\link{mboost}} to fit the components
  of an LSS (location, shape and scale) model.

  \code{blackboostLSS} uses \code{\link{blackboost}} to fit the components
  of an LSS (location, shape and scale) model, i.e., a tree-based
  boosting model is fitted for each component.

  \code{mboostLSS}, \code{glmboostLSS}, \code{gamboostLSS} and
  \code{blackboostLSS} all call \code{mboostLSS_fit} while \code{fun} is
  the corresponding \code{\link{mboost}} function, i.e., the same
  function without \code{LSS}. For further possible arguments see
  these functions as well as \code{\link{mboost_fit}}.

  In all four fitting functions it is possible to specify one or
  multiple \code{mstop} and \code{nu} values via
  \code{\link{boost_control}}. In the case of one value, this value is
  used for all parameters of the LSS model. Alternatively, a named
  list with seperate values for each component can be used to specify a
  seperate value for each parameter of the LSS model. The names of the
  list must correspond to the names of the parameters of the LSS family.

}
\value{
  An object of class \code{mboostLSS} with corresponding methods to
  extract information.
}
\references{
  Put some refs here!
}

\seealso{
  \code{\link{mboost}}, \code{\link{gamboost}}, \code{\link{glmboost}},
  \code{\link{blackboost}} for the underlying boosting functions.

  See for example \code{\link{risk}} for methods that can be used to
  extract information from \code{\link{mboostLSS}} objects.
}
\examples{
### Data generating process:
set.seed(1907)
x1 <- rnorm(1000)
x2 <- rnorm(1000)
x3 <- rnorm(1000)
x4 <- rnorm(1000)
x5 <- rnorm(1000)
x6 <- rnorm(1000)
mu    <- exp(1.5 +1 * x1 +0.5 * x2 -0.5 * x3 -1 * x4)
sigma <- exp(-0.4 * x3 -0.2 * x4 +0.2 * x5 +0.4 * x6)
y <- numeric(1000)
for( i in 1:1000)
    y[i] <- rnbinom(1, size = sigma[i], mu = mu[i])
dat <- data.frame(x1, x2, x3, x4, x5, x6, y)



### simple model with y ~ . for both components:
model <- glmboostLSS(y ~ ., families = NBinomialLSS(), data = dat,
                     control = boost_control(mstop = 400),
                     center = TRUE)
coef(model, off2int = TRUE)



### estimate model with different formulas for mu and sigma:
names(NBinomialLSS())      # names of the family

# Note: Multiple formulas must be specified via a _named list_
#       where the names correspond to the names of the family
#       (see above)
model2 <- glmboostLSS(formula = list(mu = y ~ x1 + x2 + x3 + x4,
                                    sigma = y ~ x3 + x4 + x5 + x6),
                     families = NBinomialLSS(), data = dat,
                     control = boost_control(mstop = 400, trace = TRUE),
                     center = TRUE)
coef(model2, off2int = TRUE)



### Offset needs to be specified via family:
model <- glmboostLSS(y ~ ., data = dat,
                     families = NBinomialLSS(mu = mean(mu),
                                             sigma = mean(sigma)),
                     control = boost_control(mstop = 10),
                     center = TRUE)
# Note: mu-offset = log(mean(mu)) and sigma-offset = log(mean(sigma))
#       as we use a log-link in both families
coef(model)
log(mean(mu))
log(mean(sigma))

\dontrun{## (only as is takes some time)
### use different mstop values
model3 <- glmboostLSS(formula = list(mu = y ~ x1 + x2 + x3 + x4,
                                    sigma = y ~ x3 + x4 + x5 + x6),
                     families = NBinomialLSS(), data = dat,
                     control = boost_control(mstop = list(mu = 400,
                                                          sigma = 300),
                                             trace  = TRUE),
                     center = TRUE)
coef(model3, off2int = TRUE)
}
### Alternatively we can subset model2:
# here it is assumed that the first element in the vector corresponds to
# the first parameter of model2 etc.
model2[c(400, 300)]
par(mfrow = c(1,2))
plot(model2, xlim = c(0, max(mstop(model2))))
## all.equal(coef(model2), coef(model3)) # same!
}

\keyword{models}
\keyword{nonlinear}
