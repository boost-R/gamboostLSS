\name{AftLSS}
\alias{AftLSS}

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
  GAMLSS families for accelerated failure time models
}
\description{
  Log-normal, log-logistic and Weibull accelerated failure time models for GAMLSS,
  where both distribution parameters, location (mu) and scale (sigma), are modeled via boosting.
}
\usage{
LogNormalLSS(mu = NULL, sigma = NULL)

## 'Sub-families':
LogNormalMu(mu = NULL, sigma = NULL)
LogNormalSigma(mu = NULL, sigma = NULL)


LogLogLSS(mu = NULL, sigma = NULL)

## 'Sub-families':
LogLogMu(mu = NULL, sigma = NULL)
LogLogSigma(mu = NULL, sigma = NULL)


WeibullLSS(mu = NULL, sigma = NULL)

## 'Sub-families':
WeibullMu(mu = NULL, sigma = NULL)
WeibullSigma(mu = NULL, sigma = NULL)
}

\arguments{
  \item{mu}{ offset value for distribution parameter mu . }
  \item{sigma}{ offset value for distribution parameter sigma. }
}
\details{

  Offsets can be either scalar, a vector with length equal to the number
  of observations or NULL (default). In the latter case, a scalar offset
  for this component is computed by minimizing the risk function w.r.t.
  the corresponding component (keeping the other component fixed).

}
\value{
  An object of type class \code{families}.
}

\seealso{
  Other families: \code{\link{Families}} (with constructor \code{\link{Families}})
}
\references{

Klein, J. P. and M. L. Moeschberger (2003). Survival Analysis: Techniques for Censored and
Truncated Data (Second ed.). Springer.

Mayr, A., Fenske, N., Hofner, B., Kneib, T. and Schmid, M. (2010): GAMLSS
for high-dimensional data - A flexible approach based on boosting.
Department of Statistics: Technical Reports, Nr. 98.

Rigby, R. A. and D. M. Stasinopoulos (2005). Generalized additive models for location, scale and
shape (with discussion). Applied Statistics 54, 507-554.


}
\examples{
LogNormalLSS
LogLogLSS()
WeibullLSS()
}

\keyword{ models }
