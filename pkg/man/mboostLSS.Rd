\name{mboostLSS}
\alias{mboostLSS}
\alias{blackboostLSS}
\alias{glmboostLSS}
\alias{gamboostLSS}
\alias{mboostLSS_fit}

\title{
  Fitting GAMLSS by Boosting
}
\description{
  Functions for fitting GAMLSS (generalized additive models for
  location, scale and shape) using boosting techniques. The algorithm
  circles between the distribution parameters, updating one while using
  the current fits of the others as offsets (for details see Mayr et
  al., 2010).
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
    the same formula is used for all distribution parameters. \code{formula}
    can also be a (named) list, where each list element corresponds
    to one distribution parameter of the GAMLSS distribution. The names must be 
    the same as in the family (see example for details). }
  \item{data}{ a data frame containing the variables in the model.}
  \item{families}{ an object of class \code{families}. It can be either one of
    the pre-defined distributions that come along with the package or a new distribution
    specified by the user (see \code{\link{Families}} for details).} 
  \item{control}{ a list of parameters controlling the algorithm. For
    more details see \code{\link{boost_control}}. }
  \item{weights}{ a numeric vector of weights (optional). }
  \item{fun}{ fit function. Either \code{\link{mboost}},
    \code{\link{glmboost}}, \code{\link{gamboost}} or
    \code{\link{blackboost}}. Specified directly via the according LSS
    function. E.g. \code{gamboostLSS()} calls
    \code{mboostLSS_fit(..., fun = gamboost)}. }
  \item{\dots}{Further arguments to be passed to \code{mboostLSS_fit}.
    In  \code{mboostLSS_fit}, \code{\dots} represent further arguments to be
    passed to \code{\link{mboost}} and \code{\link{mboost_fit}}. So
    \code{\dots} can be all arguments of \code{mboostLSS_fit} and
    \code{mboost_fit}. }
}
\details{

  For information on GAMLSS theory see Rigby and Stasinopoulos (2005)
  or the information provided at \url{http://gamlss.org}.

  \code{glmboostLSS} uses \code{\link[mboost]{glmboost}} to fit the
  distribution parameters of a GAMLSS -- a linear boosting model is
  fitted for each parameter.
 
  
  \code{gamboostLSS} uses \code{\link[mboost]{gamboost}} to fit the
  distribution parameters of a GAMLSS -- an additive boosting model with
  smooth effects is fitted for each parameter.

  \code{mboostLSS} uses \code{\link[mboost]{mboost}} to fit the
  distribution parameters of a GAMLSS. The type of model (linear,
  tree-based or smooth) is specified by \code{fun}.

  \code{blackboostLSS} uses \code{\link[mboost]{blackboost}} to fit the
  distribution parameters of a GAMLSS -- a tree-based boosting model is
  fitted for each parameter.

  \code{mboostLSS}, \code{glmboostLSS}, \code{gamboostLSS} and
  \code{blackboostLSS} all call \code{mboostLSS_fit} while \code{fun} is
  the corresponding \code{\link{mboost}} function, i.e., the same
  function without \code{LSS}. For further possible arguments see
  these functions as well as \code{\link{mboost_fit}}.

  In all four fitting functions it is possible to specify one or
  multiple \code{mstop} and \code{nu} values via
  \code{\link{boost_control}}. In the case of one single value, this
  value is used for all distribution parameters of the GAMLSS model.
  Alternatively, a named list with seperate values for each component
  can be used to specify a seperate value for each parameter of the
  GAMLSS model. The names of the list must correspond to the names of
  the distribution parameters of the GAMLSS family. For one-dimensional
  stopping, the user therefore can specify, e.g., \code{mstop = 100} via
  \code{\link{boost_control}}. For more-dimensional stopping, one can
  specify, e.g., \code{mstop = list(mu = 100, sigma = 200)} (see
  examples).
  

}
\value{
  An object of class \code{mboostLSS} with corresponding methods to
  extract information.
}
\references{

Mayr, A., Fenske, N., Hofner, B., Kneib, T. and Schmid, M. (2010): GAMLSS 
for high-dimensional data -- a flexible approach based on boosting. 
Department of Statistics: Technical Report 98. 
\url{http://epub.ub.uni-muenchen.de/11938/}

Rigby, R. A. and D. M. Stasinopoulos (2005). Generalized additive models for location, scale and
shape (with discussion). Applied Statistics 54, 507-554.

Buehlmann, P. and Hothorn, T. (2007), Boosting algorithms: Regularization, 
prediction and model fitting. Statistical Science, 22(4), 477--505. 
}

\seealso{
  \code{\link{mboost}}, \code{\link{gamboost}}, \code{\link{glmboost}},
  \code{\link{blackboost}} for the underlying boosting functions
  contained in \code{\link{mboost}}.

  See for example \code{\link{risk}} or \code{\link{coef}} for methods 
  that can be used to extract information from \code{\link{mboostLSS}} objects.
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



### linear model with y ~ . for both components: 400 boosting iterations
model <- glmboostLSS(y ~ ., families = NBinomialLSS(), data = dat,
                     control = boost_control(mstop = 400),
                     center = TRUE)
coef(model, off2int = TRUE)


### estimate model with different formulas for mu and sigma:
names(NBinomialLSS())      # names of the family

# Note: Multiple formulas must be specified via a _named list_
#       where the names correspond to the names of the distribution parameters
#       in the family (see above)
model2 <- glmboostLSS(formula = list(mu = y ~ x1 + x2 + x3 + x4,
                                    sigma = y ~ x3 + x4 + x5 + x6),
                     families = NBinomialLSS(), data = dat,
                     control = boost_control(mstop = 400, trace = TRUE),
                     center = TRUE)
coef(model2, off2int = TRUE)



### Offset needs to be specified via the arguments of families object:
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
### use different mstop values for the two distribution parameters
### (two-dimensional early stopping) 
### the number of iterations is passed to boost_control via a named list
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
# the first distribution parameter of model2 etc.
model2[c(400, 300)]
par(mfrow = c(1,2))
plot(model2, xlim = c(0, max(mstop(model2))))
## all.equal(coef(model2), coef(model3)) # same!

### WARNING: Subsetting via model[mstopnew] changes the model directly! 
### For the original fit one has to subset again: model[mstop]

}

\keyword{models}
\keyword{nonlinear}
\keyword{fitting}
