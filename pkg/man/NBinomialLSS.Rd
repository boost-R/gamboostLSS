\name{NBinomialLSS}
\alias{NBinomialLSS}
\alias{NBinomialMu}
\alias{NBinomialSigma}

\title{
  Family for Negative Binomial Family for LSS Models
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
NBinomialLSS(mu = 0, sigma = 1)
NBinomialMu(off = list(sigma = 1))
NBinomialSigma(off = list(mu = 0))
}

\arguments{
  \item{mu}{
    ...
}
  \item{sigma}{
    ...
}
\item{off}{
  ...
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
function (mu = 0, sigma = 1)
{
    RET <- list(mu = NBinomialMu(off = list(sigma = sigma)),
        sigma = NBinomialSigma(off = list(mu = mu)))
    class(RET) <- "families"
    return(RET)
  }
}

\keyword{ models }