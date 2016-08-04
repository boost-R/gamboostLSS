\name{weighted.median}

\alias{weighted.median}

\title{
  Weighted Median
}
\description{
  Function to compute the weighted median.
}
\usage{
weighted.median(x, w = 1, na.rm = FALSE)
}
\arguments{
  \item{x}{
    a numeric vector containing the values whose median is to be
    computed.
  }
  \item{w}{
    weights that are used to compute the median. This can be either a
    single value (which will be used as weight for all observations)
    or a numeric vector of the same length as \code{x}.
  }
  \item{na.rm}{
    logical. Should \code{NA}s (from weights \code{w} and/or data
    \code{x}) be removed?
  }
}
\details{
  The weighted median is computed as the value where the cumulative
  relative weights (relative to the sum of all weights) crosses 0.5.

  This function is used in the stabilization of the negative gradient
  via the meadian absolute deviation (MAD). For details see Hofner et al
  (2015).
}

\references{
  B. Hofner, A. Mayr, M. Schmid (2015). gamboostLSS: An R Package for
  Model Building and Variable Selection in the GAMLSS Framework.
  Journal of Statistical Software. Accepted for publication.

  Available as \code{vignette("gamboostLSS_Tutorial")}.
}

\seealso{

 \code{\link{glmboostLSS}}, \code{\link{gamboostLSS}} and
 \code{\link{blackboostLSS}} for fitting of GAMLSS where the
 standardization is explained in more detail.

}
\examples{

## compute the weighted median with case weights
x <- c(1, 2, 3, 4)
w <- c(0, 1, 2, 3)
weighted.median(x, w)

## compute the weighted median with arbitrary weights
x <- rnorm(100)
w <- runif(100)
weighted.median(x, w)

}
\keyword{methods}
