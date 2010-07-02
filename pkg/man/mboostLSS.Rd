\name{mboostLSS}
\alias{mboostLSS}
\alias{blackboostLSS}
\alias{glmboostLSS}
\alias{gamboostLSS}
\alias{mboostLSS_fit}

\title{
  Fitting mboostLSS models
}
\description{
  Functions for fitting LSS (location, shape and scale) models using
  boosting techniques.
}
\usage{
mboostLSS(...)
glmboostLSS(...)
gamboostLSS(...)
blackboostLSS(...)
mboostLSS_fit(formula, data, families = list(), ctrl = boost_control(),
              weights = NULL, fun = mboost, ...)
}

\arguments{
  \item{formula}{ a symbolic description of the model to be fit. See
    \code{\link{mboost}} for details. If formula is a single formula,
    the same formula is used for all model components.
    ADD: Can also be a (named) list, where each list element corresponds
    to one component of the LSS model. The names must be the same as in
    the family (see example for details).
}
  \item{data}{ a data frame containing the variables in the model.}
  \item{families}{ a list of mboost-\code{\link{Family}}s. }
  \item{ctrl}{ a list of parameters controlling the algorithm. For
    more details see \code{\link{boost_control}}. }
  \item{weights}{ a numeric vector of weights (optional). }
  \item{fun}{ fit function. Either \code{\link{mboost}},
    \code{\link{glmboost}}, \code{\link{gamboost}} or
    \code{\link{blackboost}}. Specified directly via the according LSS
    function. E.g. \code{gamboostLSS()} calls
    \code{mboostLSS_fit(..., fun = gamboost)}. }  
  \item{\dots}{ Further arguments to be passed to \code{\link{mboost}}
    and \code{\link{mboost_fit}}. See there for details. }
}
\details{
  Add details later.
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

\seealso{
  \code{\link{mboost}}
}
\examples{
## names of the family
names(NBinomialLSS())
}

\keyword{models}
\keyword{nonlinear}
