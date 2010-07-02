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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (mu = 0, sigma = 1, df = 1) 
{
    RET <- list(mu = StudentTMu(off = list(sigma = sigma, df = df)), 
        sigma = StudentTSigma(off = list(mu = mu, df = df)), 
        df = StudentTDf(off = list(mu = mu, sigma = sigma)))
    class(RET) <- "families"
    return(RET)
  }
}

\keyword{ models }