\name{Families}
\alias{Families}
\alias{families}

\alias{GaussianLSS}
\alias{GaussianMu}
\alias{GaussianSigma}

\alias{GammaLSS}
\alias{GammaMu}
\alias{GammaSigma}

\alias{BetaLSS}
\alias{BetaMu}
\alias{BetaPhi}


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

\alias{ZIPoLSS}
\alias{ZINBLSS}

\alias{options}
\alias{stab_ngrad}
\alias{stabilize_ngrad}
\alias{stabilize_ngradient}

\title{
  Families for GAMLSS models
}
\description{
 The package provides some pre-defined GAMLSS families, e.g.
 \code{NBionomialLSS}. Objects of the class \code{families} provide a
 convenient way to specify GAMLSS distributions to be fitted by one of
 the boosting algorithms implemented in this package. By using the
 function \code{Families}, a new object of the class \code{families}
 can be generated.
}

\usage{
############################################################
# Families for continuous response

# Gaussian distribution
GaussianLSS(mu = NULL, sigma = NULL,
            stabilization = c("none", "MAD"))

# Student's t-distribution
StudentTLSS(mu = NULL, sigma = NULL, df = NULL,
            stabilization = c("none", "MAD"))

############################################################
# Families for continuous non-negative response

# Gamma distribution
GammaLSS(mu = NULL, sigma = NULL,
         stabilization = c("none", "MAD"))

############################################################
# Families for fractions and bounded continuous response

# Beta distribution
BetaLSS(mu = NULL, phi = NULL,
        stabilization = c("none", "MAD"))

############################################################
# Families for count data

# Negative binomial distribution
NBinomialLSS(mu = NULL, sigma = NULL,
             stabilization = c("none", "MAD"))

# Zero-inflated Poisson distribution
ZIPoLSS(mu = NULL, sigma = NULL,
        stabilization = c("none", "MAD"))

# Zero-inflated negative binomial distribution
ZINBLSS(mu = NULL, sigma = NULL, nu = NULL,
        stabilization = c("none", "MAD"))

############################################################
# Families for survival models (accelerated failure time
# models) for data with right censoring

# Log-normal distribution
LogNormalLSS(mu = NULL, sigma = NULL,
             stabilization = c("none", "MAD"))

# Log-logistic distribution
LogLogLSS(mu = NULL, sigma = NULL,
          stabilization = c("none", "MAD"))

# Weibull distribution
WeibullLSS(mu = NULL, sigma = NULL,
           stabilization = c("none", "MAD"))

############################################################
# Constructor function for new GAMLSS distributions
Families(..., qfun = NULL, name = NULL)
}

\arguments{
  \item{\dots}{ sub-families to be passed to constructor. }
  \item{qfun}{ quantile function. This function can for example be used
    to compute (marginal) prediction intervals. See
    \code{\link{predint}}.}
  \item{name}{ name of the families. }
  \item{mu}{ offset value for mu. }
  \item{sigma}{ offset value for sigma. }
  \item{phi}{ offset value for phi. }
  \item{df}{ offset value for df. }
  \item{nu}{ offset value for nu. }
  \item{stabilization}{ governs if the negative gradient should be
    standardized in each boosting step. It can be either "none" or
    "MAD". See also Details below.
  }
}
\details{

  The arguments of the families are the offsets for each distribution
  parameter. Offsets can be either scalar, a vector with length equal to
  the number of observations or \code{NULL} (default). In the latter
  case, a scalar offset for this component is computed by minimizing the
  risk function w.r.t. the corresponding distribution parameter (keeping
  the other parameters fixed).

  Note that \code{gamboostLSS} is not restricted to three components but
  can handle an arbitrary number of components (which, of course,
  depends on the GAMLSS distribution). However, it is important that the
  names (for the offsets, in the sub-families etc.) are chosen
  \emph{consistently}.

  The \code{ZIPoLSS} families can be used to fit zero-inflated Poisson
  models. Here, \code{mu} and \code{sigma} refer to the location
  parameter of the Poisson component (with log link) and the mean of the
  zero-generating process (with logit link), respectively.

  Similarly, \code{ZINBLSS} can be used to fit zero-inflated negative
  binomial models. Here, \code{mu} and \code{sigma} refer to the
  location and scale parameters (with log link) of the negative binomial
  component of the model. The zero-generating process (with logit link)
  is represented by \code{nu}.

  The \code{Families} function can be used to implements a new GAMLSS
  distribution which can be used for fitting by \code{\link{mboostLSS}}.
  Thereby, the function builds a list of sub-families, one for each
  distribution parameter. The sub-families themselves are objects of the
  class \code{boost_family}, and can be constructed via the function
  \code{\link[mboost]{Family}} of the \code{mboost} Package.

  Arguments to be passed to \code{Family}: The \code{loss} for every
  distribution parameter (contained in objects of class
  \code{boost_family}) is the negative log-likelihood of the
  corresponding distribution. The \code{ngradient} is the negative
  partial derivative of the loss function with respect to the
  distribution parameter. For a two-parameter distribution (e.g. mu and
  sigma), the user therefore has to specify two sub-families with
  \code{\link{Family}}. The \code{loss} is basically the same function
  for both paramters, only \code{ngradient} differs. Both sub-families
  are passed to the \code{Families} constructor, which returns an object
  of the class \code{families}.

  To (potentially) stabilize the model estimation by standardizing the
  negative gradients one can use the argument \code{stabilization} of
  the families. If \code{stabilization = "MAD"}, the negative gradient
  is divided by its (weighted) median absolute deviation \deqn{median_i
  (|u_{k,i} - median_j(u_{k,j})|)} in each boosting step. See Hofner et
  al. (2015) for details.

}
\value{
  An object of class \code{families}.
}

\author{
  \code{BetaLSS} for boosting beta regression was implmented by Florian
  Wickler.
}

\references{
  B. Hofner, A. Mayr, M. Schmid (2015). gamboostLSS: An R Package for
  Model Building and Variable Selection in the GAMLSS Framework.
  Journal of Statistical Software. Accepted for publication.

  Available as \code{vignette("gamboostLSS_Tutorial")}.

  Mayr, A., Fenske, N., Hofner, B., Kneib, T. and Schmid, M. (2012):
  Generalized additive models for location, scale and shape for
  high-dimensional data - a flexible approach based on boosting. Journal
  of the Royal Statistical Society, Series C (Applied Statistics) 61(3):
  403-427.

  Rigby, R. A. and D. M. Stasinopoulos (2005). Generalized additive
  models for location, scale and shape (with discussion). Journal of the
  Royal Statistical Society, Series C (Applied Statistics), 54, 507-554.

}


\seealso{
  \code{\link{as.families}} for applying GAMLSS distributions provided
  in the framework of the \code{gamlss} package.

  The functions \code{\link{gamboostLSS}} and \code{\link{glmboostLSS}}
  can be used for model fitting.

  See also the corresponding constructor function
  \code{\link[mboost]{Family}} in \code{\link[mboost]{mboost}}.
}

\examples{
## Example to define a new distribution:
## Students t-distribution with two parameters, df and mu:

## sub-Family for mu
## -> generate object of the class family from the package mboost
newStudentTMu  <- function(mu, df){

    # loss is negative log-Likelihood, f is the parameter to be fitted with
    # id link -> f = mu
    loss <- function(df,  y, f) {
        -1 * (lgamma((df + 1)/2)  - lgamma(1/2) -
              lgamma(df/2) - 0.5 * log(df) -
              (df + 1)/2 * log(1 + (y - f)^2/(df )))
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

## sub-Family for df
newStudentTDf <- function(mu, df){

    # loss is negative log-Likelihood, f is the parameter to be fitted with
    # log-link: exp(f) = df
    loss <- function( mu, y, f) {
        -1 * (lgamma((exp(f) + 1)/2)  - lgamma(1/2) -
              lgamma(exp(f)/2) - 0.5 * f -
              (exp(f) + 1)/2 * log(1 + (y - mu)^2/(exp(f) )))
    }
    # risk is sum of loss
    risk <- function(y, f, w = 1) {
        sum(w * loss(y = y, f = f,  mu = mu))
    }
    # ngradient is the negative derivate w.r.t. df
    ngradient <- function(y, f, w = 1) {
        exp(f)/2 * (digamma((exp(f) + 1)/2) - digamma(exp(f)/2)) -
            0.5 - (exp(f)/2 * log(1 + (y - mu)^2 / (exp(f) )) -
                   (y - mu)^2 / (1 + (y - mu)^2 / exp(f)) * (exp(-f) + 1)/2)
    }
    # use the Family constructor of mboost
    Family(ngradient = ngradient, risk = risk, loss = loss,
           response = function(f) exp(f),
           name = "Student's t-distribution: df (log link)")
}

## families object for new distribution
newStudentT <- Families(mu= newStudentTMu(mu=mu, df=df),
                        df=newStudentTDf(mu=mu, df=df))

\donttest{### Do not test the following code per default on CRAN as it takes some time to run:
### usage of the new Student's t distribution:
library(gamlss)   ## required for rTF
set.seed(1907)
n <- 5000
x1  <- runif(n)
x2 <- runif(n)
mu <- 2 -1*x1 - 3*x2
df <- exp(1 + 0.5*x1 )
y <- rTF(n = n, mu = mu, nu = df)

## model fitting
model <- glmboostLSS(y ~ x1 + x2, families = newStudentT,
                     control = boost_control(mstop = 100),
                     center = TRUE)
## shrinked effect estimates
coef(model, off2int = TRUE)

## compare to pre-defined three parametric t-distribution:
model2 <- glmboostLSS(y ~ x1 + x2, families = StudentTLSS(),
                      control = boost_control(mstop = 100),
                      center = TRUE)
coef(model2, off2int = TRUE)

## with effect on sigma:
sigma <- 3+ 1*x2
y <- rTF(n = n, mu = mu, nu = df, sigma=sigma)
model3 <- glmboostLSS(y ~ x1 + x2, families = StudentTLSS(),
                      control = boost_control(mstop = 100),
                      center = TRUE)
coef(model3, off2int = TRUE)
}

}
\keyword{ models }
\keyword{ distributions }


