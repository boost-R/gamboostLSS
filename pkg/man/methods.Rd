\name{methods}

\alias{coef.mboostLSS}
\alias{coef.glmboostLSS}

\alias{risk}
\alias{risk.mboostLSS}

\alias{[.mboostLSS}

\alias{mstop.mboostLSS}
\alias{mstop.oobag}

\alias{selected.mboostLSS}

\alias{fitted.mboostLSS}
\alias{predict.mboostLSS}

\alias{plot.glmboostLSS}
\alias{plot.gamboostLSS}

\title{
  Methods for mboostLSS
}
\description{
  Methods for LSS models fitted by boosting algorithms.
}
\usage{
\method{coef}{mboostLSS}(object, which = NULL,
                         aggregate = c("sum", "cumsum", "none"),
                         parameter = names(object), ...)
\method{coef}{glmboostLSS}(object, which = NULL,
                           aggregate = c("sum", "cumsum", "none"),
                           off2int = FALSE, parameter = names(object), ...)

\method{risk}{mboostLSS}(object, merge = FALSE, parameter = names(object), ...)

\method{mstop}{mboostLSS}(object, parameter = names(object), ...)
\method{mstop}{oobag}(object, parameter = names(object), ...)

\method{[}{mboostLSS}(x, i, return = TRUE, ...)

\method{selected}{mboostLSS}(object)

\method{plot}{glmboostLSS}(x, main = names(x), parameter = names(x),
                           off2int = FALSE, ...)
\method{plot}{gamboostLSS}(x, main = names(x), parameter = names(x),
...)

\method{fitted}{mboostLSS}(object, parameter = names(object), ...)
\method{predict}{mboostLSS}(object, newdata = NULL,
                            type = c("link", "response", "class"),
                            which = NULL,
                            aggregate = c("sum", "cumsum", "none"),
                            parameter = names(object), ...)
}

\arguments{
  \item{x}{ an object of the appropriate class (see usage). }
  \item{object}{ an object of the appropriate class (see usage). }
  \item{which}{ a subset of base-learners to take into account for computing
     predictions or coefficients. If \code{which} is given
     (as an integer vector or characters corresponding
     to base-learners) a list or matrix is returned.}
   \item{aggregate}{ a character specifying how to aggregate predictions
     or coefficients of single base-learners. The default
     returns the prediction or coefficient for the final number of
     boosting iterations. \code{"cumsum"} returns a
     matrix with the predictions for all iterations
     simultaneously (in columns). \code{"none"} returns a
     list with matrices where the \eqn{j}th columns of the
     respective matrix contains the predictions
     of the base-learner of the \eqn{j}th boosting
     iteration (and zero if the base-learner is not
     selected in this iteration).}
  \item{parameter}{ This can be either a vector of indizes or a vector
     of parameter names which should be processed. See expamles for
     details. Per default all parameters of the LSS family are
     returned. }
  \item{off2int}{ logical indicating whether the offset should be
    added to the intercept (if there is any) or if the offset is
    neglected for plotting (default).}
  \item{merge}{ logical. Should the risk vectors of the single
    components be merged to one risk vector? Per default (\code{merge =
    FALSE}) a (named) list of risk vectors is returned. }
  \item{i}{ integer. Index specifying the model to extract. If \code{i}
    is smaller than the initial \code{mstop}, a subset is used.
    If \code{i} is larger than the initial \code{mstop},
    additional boosting steps are performed until step \code{i}
    is reached. See details for more information. }
  \item{return}{ a logical indicating whether the changed object is
    returned. }
  \item{main}{ a title for the plots.}
   \item{newdata}{ optionally; A data frame in which to look for variables with
    which to predict.}
  \item{type}{ the type of prediction required.  The default is on the scale
    of the predictors; the alternative \code{"response"} is on
    the scale of the response variable.  Thus for a
    binomial model the default predictions are of log-odds
    (probabilities on logit scale) and \code{type = "response"} gives
    the predicted probabilities.  The \code{"class"} option returns
    predicted classes.}
  \item{\dots}{ Further arguments to the functions. }
}
\details{
  These functions can be used to extract details from fitted models.
  \code{print} shows a dense representation of the model fit.

  ...
}

\section{Warning}{
  The \code{[.mboostLSS} function changes the original object, i.e.
  \code{LSSmodel[10]} changes \code{LSSmodel} directly!
}

\seealso{
  ...
}
\examples{

### generate data
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

### fit a model
model <- gamboostLSS(y ~ ., families = NBinomialLSS(), data = dat,
                     control = boost_control(mstop = 400))

### extract coefficients
coef(model)
### plot complete model
par(mfrow=c(4,3))
plot(model)
### plot first parameter only
par(mfrow=c(2,3))
plot(model, parameter = "mu")
### now plot only effect of x3 of both parameters
par(mfrow = c(1, 2))
plot(model, which = "x3")
### first component second parameter (sigma)
par(mfrow = c(1, 1))
plot(model, which = 1, parameter = 2)

}
\keyword{methods}

