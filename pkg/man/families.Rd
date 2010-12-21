\name{Families}
\alias{Families}

\alias{NBinomialLSS}
\alias{NBinomialMu}
\alias{NBinomialSigma}

\alias{StudentTLSS}
\alias{StudentTMu}
\alias{StudentTSigma}
\alias{StudentTDf}

\alias{LogNormalLSS}
\alias{LogNormalMu}
\alias{LogNormalSigma}

\alias{WeibullLSS}
\alias{WeibullMu}
\alias{WeibullSigma}

\alias{LogLogLSS}
\alias{LogLogMu}
\alias{LogLogSigma}


\title{
  Families for GAMLSS models
}
\description{

 Objects of the class \code{families} provide a convenient way to specify GAMLSS 
 distributions to be fitted by one of the boosting 
 algorithms implemented in this package. The package provides some pre-defined GAMLSS families, 
 but the user can also specify new distributions using the \code{Families} function. 
 
}
\usage{
## Constructor function for new GAMLSS distributions
Families(...)

###  Pre-defined GAMLSS distributions  ###

## families Object for the negative binomial distribution 
NBinomialLSS(mu = NULL, sigma = NULL)
## 'Sub-families' for the two distribution parameters:
NBinomialMu(mu = NULL, sigma = NULL)
NBinomialSigma(mu = NULL, sigma = NULL)

## families Object for Student's t-distribution 
StudentTLSS(mu = NULL, sigma = NULL, df = NULL)
## 'Sub-families' for the three distribution parameters:
StudentTMu(mu = NULL, sigma = NULL, df = NULL)
StudentTSigma(mu = NULL, sigma = NULL, df = NULL)
StudentTDf(mu = NULL, sigma = NULL, df = NULL)

## families Object for log-normal accelerated failure time model
LogNormalLSS(mu = NULL, sigma = NULL)
## 'Sub-families' for the two distribution parameters:
LogNormalMu(mu = NULL, sigma = NULL)
LogNormalSigma(mu = NULL, sigma = NULL)

## families Object for log-logistic accelerated failure time model
LogLogLSS(mu = NULL, sigma = NULL)
## 'Sub-families' for the two distribution parameters:
LogLogMu(mu = NULL, sigma = NULL)
LogLogSigma(mu = NULL, sigma = NULL)

## families Object for accelerated failure time model with Weibull distribution
WeibullLSS(mu = NULL, sigma = NULL)
## 'Sub-families' for the two distribution parameters:
WeibullMu(mu = NULL, sigma = NULL)
WeibullSigma(mu = NULL, sigma = NULL)
}

\arguments{
  \item{mu}{ offset value for mu. }
  \item{sigma}{ offset value for sigma. }
  \item{df}{ offset value for df. }
  \item{\dots}{ sub-families to be passed to constructor. }
}
\details{
  
  \code{Families} implements a new GAMLSS distribution to be used for
  fitting by \code{\link{mboostLSS}}. The function builds a list of
  Sub-families, one for each distribution parameter. The Sub-families
  itself are objects of the class \code{boost_family}, and can be
  constructed via the function \code{\link[mboost]{Family}} of the
  \code{mboost} Package. The \code{loss} for every distribution
  parameter (contained in objects of class \code{boost_family}) is the
  negative log-likelihood of the corresponding distribution. The
  \code{ngradient} is the negative partial derivative
  of the loss function with respect to the distribution parameter. For a
  two-parameter distribution (e.g. mu and sigma), the user therefore has
  to specify two sub-families with \code{\link{Family}}. The \code{loss}
  is basically the same function for both paramters, only the
  \code{ngradient} differs. Both sub-families are passed to the
  \code{Families} constructor, which returns in a \code{families}
  object.
    
  The arguments of the final object are the offsets for each
  distribution parameter. Offsets can be either scalar, a vector with
  length equal to the number of observations or NULL (default). In the
  latter case, a scalar offset for this component is computed by
  minimizing the risk function w.r.t. the corresponding distribution
  parameter (keeping the other parameters fixed).

  Note that gamboostLSS is not restricted to the three components
  \code{mu}, \code{sigma} and \code{df} but can handle an arbitrary
  number of components (which, of course, depends on the GAMLSS distribution).
  However, it is important that the names (for the offsets, in the subfamilies
  etc.) are chosen \emph{consistently}. }
  \value{ An object of class \code{families}. }

\references{
Mayr, A., Fenske, N., Hofner, B., Kneib, T. and Schmid, M. (2010): GAMLSS 
for high-dimensional data - a flexible approach based on boosting. 
Department of Statistics: Technical Report 98. \url{http://epub.ub.uni-muenchen.de/11938/}

Rigby, R. A. and D. M. Stasinopoulos (2005). Generalized additive models for location, scale and
shape (with discussion). Applied Statistics 54, 507-554.
}

\seealso{
 
 The corresponding constructor function \code{\link[mboost]{Family}} in
 \code{\link[mboost]{mboost}}. For the fitting of GAMLSS distributions
 see \code{\link{mboostLSS}}.
}
\examples{

# Example to define a new distribution:
# Students t-distribution with two parameters, df and mu:

# Sub-Family for mu
newStudentTMu  <- function(mu, df){

    # loss is negative log-Likelihood, f is the parameter to be fitted with 
    # id link -> f = mu
    loss <- function(df,  y, f) {
        -1 * (lgamma((df + 1)/2)  - lgamma(1/2) -
            lgamma(df/2) - 0.5 * log(df) - (df + 1)/2 * log(1 +
            (y - f)^2/(df )))
    }
    # risk is sum of loss
    risk <- function(y, f, w = 1) {
        sum(w * loss(y = y, f = f, df = df))
    }
    # ngradient is the negative derivate w.r.t. mu
    ngradient <- function(y, f, w = 1) {
        (df + 1) * (y - f)/(df  + (y - f)^2)
    }

    # use the Family constructor of mboost
    Family(ngradient = ngradient, risk = risk, loss = loss,
        response = function(f) f,  
        name = "new Student's t-distribution: mu (id link)")
}

# Sub-Family for df
newStudentTDf <- function(mu, df){

     # loss is negative log-Likelihood, f is the parameter to be fitted with
     # log-link: exp(f) = df
    loss <- function( mu, y, f) {
        -1 * (lgamma((exp(f) + 1)/2)  - lgamma(1/2) -
            lgamma(exp(f)/2) - 0.5 * f - (exp(f) + 1)/2 * log(1 +
            (y - mu)^2/(exp(f) )))
    }
    # risk is sum of loss
    risk <- function(y, f, w = 1) {
        sum(w * loss(y = y, f = f,  mu = mu))
    }
        # ngradient is the negative derivate w.r.t. df
     ngradient <- function(y, f, w = 1) {
        exp(f)/2 * (digamma((exp(f) + 1)/2) - digamma(exp(f)/2)) -
            0.5 - (exp(f)/2 * log(1 + (y - mu)^2/(exp(f) )) -
            (y - mu)^2/(1 + (y - mu)^2/exp(f)) * (exp(-f) +
                1)/2)
    }
    # use the Family constructor of mboost
      Family(ngradient = ngradient, risk = risk, loss = loss,
        response = function(f) exp(f),
        name = "Student's t-distribution: df (log link)")
}
    
# families object for new distribution
newStudentT <- Families(mu= newStudentTMu(mu=mu, df=df), 
                                    df=newStudentTDf(mu=mu, df=df))


### usage of the new Student's t distribution:
library(gamlss)   ## required for rTF
set.seed(1907)
n <- 5000
x1  <- runif(n)
x2 <- runif(n)
mu <- 2 -1*x1 - 3*x2
df <- exp(1 + 0.5*x1 )
y <- rTF(n = n, mu = mu, nu = df)

# model fitting
model <- glmboostLSS(y ~ x1 + x2, families = newStudentT,
                     control = boost_control(mstop = 100),
                     center = TRUE)
# shrinked effect estimates                     
coef(model, off2int = TRUE)

# compare to pre-defined three parametric t-distribution:
model2 <- glmboostLSS(y ~ x1 + x2, families = StudentTLSS(),
                     control = boost_control(mstop = 100),
                     center = TRUE)
coef(model2, off2int = TRUE)

# with effect on sigma:
sigma <- 3+ 1*x2
y <- rTF(n = n, mu = mu, nu = df, sigma=sigma)
model3 <- glmboostLSS(y ~ x1 + x2, families = StudentTLSS(),
                     control = boost_control(mstop = 100),
                     center = TRUE)
coef(model3, off2int = TRUE)



}
\keyword{ models }
\keyword{ distributions }

