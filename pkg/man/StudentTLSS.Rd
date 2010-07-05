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
StudentTLSS(mu = 0, sigma = 1, df = 1)
StudentTMu(off = list(df = 1, sigma = 1))
StudentTSigma(off = list(df = 1, mu = 0))
StudentTDf(off = list(sigma = 1, mu = 0))


}

\arguments{
  \item{mu}{
%%     ~~Describe \code{mu} here~~
}
  \item{sigma}{
%%     ~~Describe \code{sigma} here~~
}
  \item{df}{
%%     ~~Describe \code{df} here~~
}
\item{off}{
%%     ~~Describe \code{off} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
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