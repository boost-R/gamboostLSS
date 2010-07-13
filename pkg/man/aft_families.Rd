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
  Families for AFT Models
}
\description{
  Log-normal, log-logistic and Weibull AFT models where both location
  and scale are modeled via boosting.
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
  \item{mu}{ offset value for mu component. }
  \item{sigma}{ offset value for sigma component. }
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
  Other families: \code{\link{Families}} (with contructor \code{\link{makeFamilies}})
}
\examples{
LogNormalLSS
LogLogLSS()
WeibullLSS()
}

\keyword{ models }