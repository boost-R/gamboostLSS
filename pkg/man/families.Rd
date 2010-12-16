\name{Families}
\alias{Families}

\alias{NBinomialLSS}
\alias{NBinomialMu}
\alias{NBinomialSigma}

\alias{StudentTLSS}
\alias{StudentTMu}
\alias{StudentTSigma}
\alias{StudentTDf}

\title{
  Families for GAMLSS models
}
\description{
  Negative binomial and Student's t distribution for GAMLSS models fitted
  via boosting.
}
\usage{
## Constructor function
Families(...)


NBinomialLSS(mu = NULL, sigma = NULL)

## 'Sub-families':
NBinomialMu(mu = NULL, sigma = NULL)
NBinomialSigma(mu = NULL, sigma = NULL)


StudentTLSS(mu = NULL, sigma = NULL, df = NULL)

## 'Sub-families':
StudentTMu(mu = NULL, sigma = NULL, df = NULL)
StudentTSigma(mu = NULL, sigma = NULL, df = NULL)
StudentTDf(mu = NULL, sigma = NULL, df = NULL)
}

\arguments{
  \item{mu}{ offset value for mu component. }
  \item{sigma}{ offset value for sigma component. }
  \item{df}{ offset value for df component. }
  \item{\dots}{ sub-families to be passed to constructor. }
}
\details{

  Offsets can be either scalar, a vector with length equal to the number
  of observations or NULL (default). In the latter case, a scalar offset
  for this component is computed by minimizing the risk function w.r.t.
  the corresponding component (keeping the other component fixed).

}
\value{
  An object of class \code{families}.
}

\seealso{
  AFT families for survival models: \code{\link{AftLSS}}
}
\examples{
NBinomialLSS()

### Student's t:
library(gamlss)   ## required for rTF
set.seed(1907)
n <- 5000
x1  <- runif(n)
x2 <- runif(n)
mu <- 2 -1*x1 - 3*x2
sigma <- exp(-1*x1 + 3*x2)
df <- exp(1 + 3*x1 + 1*x2)
y <- rTF(n = n, mu = mu, sigma = sigma, nu = df)

model <- glmboostLSS(y ~ x1 + x2, families = StudentTLSS(),
                     control = boost_control(mstop = 100),
                     center = TRUE)
coef(model, off2int = TRUE)
selected(model)
}

\keyword{ models }