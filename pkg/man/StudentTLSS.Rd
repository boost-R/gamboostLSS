\name{StudentTLSS}
\alias{StudentTLSS}
\alias{StudentTMu}
\alias{StudentTSigma}
\alias{StudentTDf}

\title{
  Family for Student's t Distribution for LSS Models
}
\description{
  ...
}
\usage{
StudentTLSS(mu = 0, sigma = 1, df = 4)

## 'Sub-families':
StudentTMu(mu = NULL, sigma = NULL, df = NULL)
StudentTSigma(mu = NULL, sigma = NULL, df = NULL)
StudentTDf(mu = NULL, sigma = NULL, df = NULL)
}

\arguments{
  \item{mu}{
    ...
  }
  \item{sigma}{
    ...
  }
  \item{df}{
    ...
  }
}
\details{

}
\value{
  An object of type class \code{families}.
}
\references{

}

\seealso{
  Other families: \code{\link{NBinomialLSS}}
}
\examples{
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
                     ctrl = boost_control(mstop = 100),
                     center = TRUE)
coef(model, off2int = TRUE)
selected(model)
}

\keyword{ models }