\name{as.families}
\alias{as.families}

\alias{gamlss.Families}
\alias{gamlss1parMu}

\alias{gamlss2parMu}
\alias{gamlss2parSigma}

\alias{gamlss3parMu}
\alias{gamlss3parSigma}
\alias{gamlss3parNu}

\alias{gamlss4parMu}
\alias{gamlss4parSigma}
\alias{gamlss4parNu}
\alias{gamlss4parTau}

\title{
  Include \code{gamlss} families in the boosting framework of \code{gamboostLSS}
}

\description{
  The function \code{as.families()} provides an interface to apply
  the available distributions (families) of the \code{gamlss.dist} package
  for boosting GAMLSS via \code{gamboostLSS}.

  The function automatically builds sub-families for every distribution
  parameter and uses the constructor function \code{\link{Families}} to
  build a \code{families} object, which can be then included in the
  fitting functions \code{\link{gamboostLSS}} and
  \code{\link{glmboostLSS}}.
}
\usage{
as.families(fname = "NO", stabilization = c("none", "MAD", "L2"),
            mu = NULL, sigma = NULL, nu = NULL, tau = NULL,
            mu.link = NULL, sigma.link = NULL, nu.link = NULL, 
            tau.link = NULL)

## a wrapper to as.families:
gamlss.Families(...)
}

\arguments{
  \item{fname}{
    name of the distribution in the \code{gamlss} framework,
    as specified in the \code{gamlss.dist} package (e.g., \code{"NO"} for a
    normal distribution with parameters \code{mu} and \code{sigma}).
    Alternatively, one can directly specify the function (i.e.,
    \code{NO}) or the evaluated function \code{NO()}.
  }
  \item{mu}{possible offset value for parameter \code{mu}.}
  \item{sigma}{possible offset value for parameter \code{sigma}.}
  \item{nu}{possible offset value for parameter \code{nu}.}
  \item{tau}{possible offset value for parameter \code{tau}.}
  \item{mu.link}{different link function for parameter \code{mu}.}
  \item{sigma.link}{different link function for parameter \code{sigma}.}
  \item{nu.link}{different link function for parameter \code{nu}.}
  \item{tau.link}{different link function for parameter \code{tau}.}
  \item{stabilization}{ governs if the negative gradient should be
    standardized in each boosting step. It can be either "none" or
    "MAD". For details see \code{\link{Families}}.
  }
  \item{\ldots}{same arguments as above.}
}
\details{

  The function aims at providing an interface to include all available
  GAMLSS distributions which are implemented with the original
  \code{gamlss.dist} package in the model-based boosting framework. The user
  specifies the name of the family (as it is called in
  \code{gamlss.dist}), and the function automatically builds the
  corresponding \code{mboost}-like sub-families and the final
  \code{families} object, which can be then used with the fitting
  functions \code{\link{gamboostLSS}} and \code{\link{glmboostLSS}}.
  
  If no different link functions are specified, the standard links for the
  corresponding family in \code{gamlss.dist} are applied. 

  To extract the necessary information regarding partial derivatives
  (for the \code{ngradient} - see \code{\link{Family}} for details) and
  the log-likelihood (for the \code{loss}) the \code{gamlss.dist}
  package is loaded. If the package is not installed yet, this will
  prompt an error message.

  The functions \code{gamlss1parMu}, \code{gamlss2parMu},
  \code{gamlss2parSigma}, ... , \code{gamlss4parTau} are called
  internally to construct the sub-families. For one-parametric
  distributions, the function will prompt a warning and returns a
  \code{\link{mboost}} family, which can be then used by the fitting
  functions of the \code{\link{mboost}} package.

  For information on GAMLSS theory see Rigby and Stasinopoulos (2005),
  lists of available distributions are provided at
  \url{https://www.gamlss.com/}. For more on details boosting GAMLSS see Mayr
  et al. (2012). Hofner et al. (2016) provides a worked example and more
  details on \code{as.families}.

  To (potentially) stabilize the model estimation by standardizing the
  negative gradients one can use the argument \code{stabilization} of
  the families. See \code{\link{Families}} for details.
}

\author{
  The help of Mikis Stasinpoulos during the work on this function is
  gratefully acknowledged.
}

\value{
  An object of class \code{families}. If the user specifies a
  one-parametric distribution, an object of class \code{family} is
  returned.
}

\references{

  B. Hofner, A. Mayr, M. Schmid (2016). gamboostLSS: An R Package for
  Model Building and Variable Selection in the GAMLSS Framework.
  Journal of Statistical Software, 74(1), 1-31.

  Available as \code{vignette("gamboostLSS_Tutorial")}.

  Mayr, A., Fenske, N., Hofner, B., Kneib, T. and Schmid, M. (2012):
  Generalized additive models for location, scale and shape for
  high-dimensional data - a flexible approach based on boosting. Journal
  of the Royal Statistical Society, Series C (Applied Statistics) 61(3):
  403-427.

  Rigby, R. A. and D. M. Stasinopoulos (2005). Generalized additive models
  for location, scale and shape (with discussion). Journal of the Royal
  Statistical Society, Series C (Applied Statistics), 54, 507-554.

}

\seealso{
  \code{gamlss.dist} for available distributions in the \code{gamlss}
  framework.

  \code{\link{Families}} for a documentation of pre-implemented
  distributions for \code{gamboostLSS}, as well as possibilities for
  user-defined distributions.
}
\examples{
## simulate small example
set.seed(123)
x <- runif(1000)

y <- rnorm(mean = 2 + 3 * x,        # effect on mu
           sd   = exp( 1 - 1 * x ), # effect on sigma
           n    = 1000)

## boosting
glmss <- glmboostLSS(y ~ x, families = as.families("NO"))
## the same:
if (require("gamlss.dist")) {
    glmss <- glmboostLSS(y ~ x, families = as.families(NO))
    glmss <- glmboostLSS(y ~ x, families = as.families(NO()))
}

coef(glmss, off2int = TRUE)

## compare to gamlss
library(gamlss)
glmss2 <- gamlss(y ~ x, sigma.formula = ~x, family = "NO")
coef(glmss2)
glmss2$sigma.coef

}

\keyword{models}
\keyword{distributions}

