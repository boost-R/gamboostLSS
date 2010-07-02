\name{methods}

\alias{coef.mboostLSS}

\alias{risk}
\alias{risk.mboostLSS}

\alias{[.mboostLSS}

\alias{mstop.mboostLSS}
\alias{mstop.oobag}

\title{
  Methods for mboostLSS
}
\description{
  Methods for LSS models fitted by boosting algorithms.
}
\usage{
\method{coef}{mboostLSS}(object, ...)

\method{risk}{mboostLSS}(object, merge = FALSE, ...)

\method{mstop}{mboostLSS}(object, ...)
\method{mstop}{oobag}(object, ...)

\method{[}{mboostLSS}(x, i, return = TRUE, ...)
}

\arguments{
  \item{x}{ model of class \code{mboostLSS}. }
  \item{object}{ model of class \code{mboostLSS}. }
  \item{merge}{ logical. Should the risk vectors of the single
    components be merged to one risk vector? Per default (\code{merge =
      FALSE}) a (named) list of risk vectors is returned. }
  \item{i}{ integer. Index specifying the model to extract. If \code{i}
            is smaller than the initial \code{mstop}, a subset is used.
            If \code{i} is larger than the initial \code{mstop},
            additional boosting steps are performed until step \code{i}
            is reached. See details for more information. }
  \item{return}{ a logical indicating whether the changed object is
                 returned. }
  \item{\dots}{ Further arguments to the functions. }
}
\details{
  ADD SOMETHING
}

\section{Warning}{
  The \code{[.mboostLSS} function changes the original object, i.e.
  \code{LSSmodel[10]} changes \code{LSSmodel} directly!
}


\references{
  ...
}

\seealso{
  ...
}
\examples{
### ADD SOMETHING
}
\keyword{methods}

