\name{methods}

\alias{coef.mboostLSS}
\alias{coef.glmboostLSS}

\alias{risk}
\alias{risk.mboostLSS}

\alias{[.mboostLSS}

\alias{mstop.mboostLSS}
\alias{mstop.oobag}
\alias{mstop.cvriskLSS}

\alias{selected}
\alias{selected.mboostLSS}

\alias{fitted.mboostLSS}
\alias{predict.mboostLSS}

\alias{predint}
\alias{PI}

\alias{plot.glmboostLSS}
\alias{plot.gamboostLSS}
\alias{plot.predint}

\alias{update.mboostLSS}

\alias{model.weights}
\alias{model.weights.default}
\alias{model.weights.mboostLSS}

\title{
  Methods for mboostLSS
}
\description{
  Methods for GAMLSS models fitted by boosting algorithms.
}
\usage{
## extract coefficients
\method{coef}{glmboostLSS}(object, which = NULL,
     aggregate = c("sum", "cumsum", "none"),
     off2int = FALSE, parameter = names(object), ...)
\method{coef}{mboostLSS}(object, which = NULL,
     aggregate = c("sum", "cumsum", "none"),
     parameter = names(object), ...)

## plot partial effects
\method{plot}{glmboostLSS}(x, main = names(x), parameter = names(x),
     off2int = FALSE, ...)
\method{plot}{gamboostLSS}(x, main = names(x), parameter = names(x), ...)

## extract and plot marginal prediction intervals
predint(x, which, pi = 0.9, newdata = NULL, ...)
PI(x, which, pi = 0.9, newdata = NULL, ...)
\method{plot}{predint}(x, main = "Marginal Prediction Interval(s)",
     xlab = NULL, ylab = NULL, lty = c("solid", "dashed"),
     lcol = c("black", "black"), log = "", ...)

## extract mstop
\method{mstop}{mboostLSS}(object, parameter = names(object), ...)
\method{mstop}{oobag}(object, parameter = names(object), ...)
\method{mstop}{cvriskLSS}(object, parameter = NULL, ...)

## set mstop
\method{[}{mboostLSS}(x, i, return = TRUE, ...)

## extract risk
\method{risk}{mboostLSS}(object, merge = FALSE, parameter = names(object), ...)

## extract selected base-learners
\method{selected}{mboostLSS}(object, parameter = names(object), ...)

## extract fitted values
\method{fitted}{mboostLSS}(object, parameter = names(object), ...)

## make predictions
\method{predict}{mboostLSS}(object, newdata = NULL,
        type = c("link", "response", "class"), which = NULL,
        aggregate = c("sum", "cumsum", "none"),
        parameter = names(object), ...)

## update weights of the fitted model
\method{update}{mboostLSS}(object, weights, oobweights = NULL,
       risk = NULL, trace = NULL, mstop = NULL, ...)

## extract model weights
\method{model.weights}{mboostLSS}(x, ...)

}
\arguments{
  \item{x}{ an object of the appropriate class (see usage). }
  \item{object}{ an object of the appropriate class (see usage). }
  \item{which}{ a subset of base-learners to take into account when computing
     predictions or coefficients. If \code{which} is given
     (as an integer vector or characters corresponding
     to base-learners), a list or matrix is returned. In \code{plot_PI}
     the argument \code{which} must be specified and it must be given as
     a character string containing the name of the variable.}
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
  \item{parameter}{ This can be either a vector of indices or a vector
     of parameter names which should be processed. See expamles for
     details. Per default all distribution parameters of the GAMLSS family are
     returned. }
  \item{off2int}{ logical indicating whether the offset should be
    added to the intercept (if there is any) or if the offset is
    neglected for plotting (default).}
  \item{merge}{ logical. Should the risk vectors of the single
    components be merged to one risk vector for the model in total? Per
    default (\code{merge = FALSE}) a (named) list of risk vectors is
    returned.}
  \item{i}{ integer. Index specifying the model to extract. If \code{i}
    is smaller than the initial \code{mstop}, a subset is used.
    If \code{i} is larger than the initial \code{mstop},
    additional boosting steps are performed until step \code{i}
    is reached. One can specify a scalar, a (possibly named) vector or a
    (possibly named) list with separate values  for each component. See
    the details section of \code{\link{mboostLSS}} for more information.}
  \item{return}{ a logical indicating whether the changed object is
    returned. }
  \item{main}{ a title for the plots.}
  \item{xlab, ylab}{ x- and y axis labels for the plots.}
  \item{pi}{ the level(s) of the prediction interval(s); Per default a
    90\% prediction interval is used.}
  \item{lty}{ (vector) of line types to be used for plotting the
    prediction intervals. The vector should contain \code{length(pi) +
    1} elements. If less elements are specified, the last element is
    recycled. The first value \code{lty[1]} is used for the marginal
    median, the second value \code{lty[2]} is used for the \code{pi[1]}
    prediction interval, etc.}
  \item{lcol}{ (vector) of (line) colors to be used for plotting the
    prediction intervals. The vector should contain \code{length(pi) +
    1} elements. If less elements are specified, the last element is
    recycled. The first value \code{lcol[1]} is used for the marginal
    median, the second value \code{lcol[2]} is used for the \code{pi[1]}
    prediction interval, etc.}
  \item{log}{ a character string which determines if and if so which
    axis should be logarithmic. See \code{\link{plot.default}} for details.}
  \item{newdata}{ optional; A data frame in which to look for variables with
    which to predict or with which to plot the marginal prediction intervals.}
  \item{type}{ the type of prediction required.  The default is on the scale
    of the predictors; the alternative \code{"response"} is on
    the scale of the response variable.  Thus for a
    binomial model the default predictions are on the log-odds scale
    (probabilities on logit scale) and \code{type = "response"} gives
    the predicted probabilities.  The \code{"class"} option returns
    predicted classes.}
  \item{weights}{
    a numeric vector of weights for the model
  }
  \item{oobweights}{
    an additional vector of out-of-bag weights (used internally
    by \code{\link[gamboostLSS]{cvrisk}}. For details see there.).
  }
  \item{risk}{
    a character indicating how the empirical risk should be
    computed for each boosting iteration. Per default \code{risk} is set
    to the risk type specified for model fitting via
    \code{\link{boost_control}}. For details and alternatives see there.
  }
  \item{trace}{
    a logical triggering printout of status information during the
    fitting process.
  }
  \item{mstop}{
    number of boosting iterations.
  }
  \item{\dots}{
    Further arguments to the functions.
  }
}
\details{
  These functions can be used to extract details from fitted models.
  \code{print} shows a dense representation of the model fit.

  The function \code{coef} extracts the regression coefficients of
  linear predictors fitted using the \code{\link{glmboostLSS}} function or
  additive predictors fitted using \code{\link{gamboostLSS}}. Per default,
  only coefficients of selected base-learners are returned for all
  distribution parameters. However, any desired coefficient can be
  extracted using the \code{which} argument. Furhtermore, one can
  extract only coefficients for a single distribution parameter via the
  \code{parameter} argument (see examples for details).

  Analogical, the function \code{plot} per default displays the
  coefficient paths for the complete GAMLSS but can be restricted to
  single distribution parameters or covariates (or subsets) using the
  \code{parameter} or \code{which} arguments, respectively.

  The function \code{predint} (or \code{PI} which is just an alias)
  computes marginal prediction intervals and returns a data frame with
  the predictors used for the marginal prediction interval, the computed
  median prediction and the marginal prediction intervals. A plot
  function (\code{plot.predint}) for the resulting object exists. Note
  that marginal predictions from AFT models (i.e., families
  \code{\link{LogLogLSS}}, \code{\link{LogNormalLSS}}, and
  \code{\link{WeibullLSS}}) represent the predicted \dQuote{true}
  survival time and not the observed survival time which is possible
  subject to censoring. Hence, comparing observed survival times with
  the marginal prediction interval is only sensible for uncensored
  observations.

  The \code{predict} function can be used for predictions for the
  distribution parameters depending on new observations whereas
  \code{fitted} extracts the regression fits for the observations in the
  learning sample. For \code{predict}, \code{newdata} can be specified
  -- otherwise the fitted values are returned. If \code{which} is
  specified, marginal effects of the corresponding base-learner(s) are
  returned. The argument \code{type} can be used to make predictions on
  the scale of the link (i.e., the linear predictor X * beta), the
  \code{response} (i.e. h(X * beta), where h is the response function)
  or the \code{class} (in case of classification).

  The function \code{update} updates models fit with \pkg{gamboostLSS}
  and is primarily used within \code{\link[gamboostLSS]{cvrisk}}. It
  updates the weights and refits the model to the altered data.
  Furthermore, the type of \code{risk}, the \code{trace} and the number
  of boosting iterations \code{mstop} can be modified.

  The function \code{model.weights} is a generic version of the same
  function provided by package \pkg{stats}, which is required to make
  \code{model.weights} work with \code{mboostLSS} models.

}

\section{Warning}{
  The \code{[.mboostLSS} function changes the original object, i.e.,
  \code{LSSmodel[10]} changes \code{LSSmodel} directly!
}

\references{

Mayr, A., Fenske, N., Hofner, B., Kneib, T. and Schmid, M. (2012):
Generalized additive models for location, scale and shape for
high-dimensional data - a flexible approach based on boosting. Journal
of the Royal Statistical Society, Series C (Applied Statistics) 61(3):
403-427.

Buehlmann, P. and Hothorn, T. (2007), Boosting algorithms:
regularization, prediction and model fitting. Statistical Science,
22(4), 477--505.

Rigby, R. A. and D. M. Stasinopoulos (2005). Generalized additive models
for location, scale and shape (with discussion). Journal of the Royal
Statistical Society, Series C (Applied Statistics), 54, 507-554.

}

\seealso{

 \code{\link{glmboostLSS}}, \code{\link{gamboostLSS}} and
 \code{\link{blackboostLSS}} for fitting of GAMLSS.

 Available distributions (families) are documented here:
 \code{\link{Families}}.

 See \code{\link[mboost]{methods}} in the \code{mboost} package for the
 corresponding methods for \code{\link[mboost]{mboost}} objects.

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
                     control = boost_control(mstop = 100))

\donttest{### use a model with more iterations for a better fit
model[400]
}
### extract coefficients
coef(model)

### only for distribution parameter mu
coef(model, parameter = "mu")

### only for covariate x1
coef(model, which = "x1")


### plot complete model
par(mfrow = c(4, 3))
plot(model)
### plot first parameter only
par(mfrow = c(2, 3))
plot(model, parameter = "mu")
### now plot only effect of x3 of both parameters
par(mfrow = c(1, 2))
plot(model, which = "x3")
### first component second parameter (sigma)
par(mfrow = c(1, 1))
plot(model, which = 1, parameter = 2)

\donttest{### plot marginal prediction interval
pi <- predint(model, pi = 0.9, which = "x1")
pi <- predint(model, pi = c(0.8, 0.9), which = "x1")
plot(pi, log = "y")  # warning as some y values are below 0
## here it would be better to plot x1 against
## sqrt(y) and sqrt(pi)

### subset model for mstop = 300 (one-dimensional)
model[300]
}
# WARNING: Subsetting via model[mstopnew] changes the model directly!
# For the original fit one has to subset again: model[mstop]

par(mfrow = c(2, 2))
plot(risk(model, parameter = "mu")[[1]])
plot(risk(model, parameter = "sigma")[[1]])

\donttest{### get back to orignal fit
model[400]
plot(risk(model, parameter = "mu")[[1]])
plot(risk(model, parameter = "sigma")[[1]])

### use different mstop values for the components
model[c(100, 200)]
## same as
  model[c(mu = 100, sigma = 200)]
## or
  model[list(mu = 100, sigma = 200)]
## or
  model[list(100, 200)]

plot(risk(model, parameter = "mu")[[1]])
plot(risk(model, parameter = "sigma")[[1]])
}
}
\keyword{methods}
