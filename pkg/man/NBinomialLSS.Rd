\name{NBinomialLSS}
\alias{NBinomialLSS}
\alias{NBinomialMu}
\alias{NBinomialSigma}

\title{
  Family for Negative Binomial Family for LSS Models
}
\description{
  ...
}
\usage{
NBinomialLSS(mu = NULL, sigma = NULL)

## 'Sub-families':
NBinomialMu(mu = NULL, sigma = NULL)
NBinomialSigma(mu = NULL, sigma = NULL)
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
  Other families: \code{\link{StudentTLSS}}
}
\examples{
NBinomialLSS()
}

\keyword{ models }