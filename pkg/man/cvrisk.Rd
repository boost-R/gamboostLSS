\name{cvrisk}

\alias{cvrisk.mboostLSS}
\alias{make.grid}

\title{ Cross-Validation }
\description{
  Cross-validated estimation of the empirical risk for hyper-parameter selection.
}

\usage{
\method{cvrisk}{mboostLSS}(object, folds = cv(model.weights(object)),
       grid = make.grid(mstop(object)),
       papply = mclapply,
       fun = NULL, ...)

make.grid(max, length.out = 10, min = NULL, log = TRUE)


}
\arguments{
  \item{object}{
    an object of class \code{mboostLSS}, i.e. a boosted GAMLSS model.
  }
  \item{folds}{
    a weight matrix with number of rows equal to the number of
    observations. The number of columns corresponds to the number of
    cross-validation runs. Can be computed using function \code{cv} and
    defaults to 25 bootstrap samples.
  }
  \item{grid}{
    a list of stopping parameters the empirical risk is to be evaluated
    for. See \code{make.grid} for an easy way to construct the grid.
  }
  \item{papply}{
    (parallel) apply function, defaults to
    \code{\link[parallel]{mclapply}}. Alternatively, \code{parLapply}
    can be used. In the latter case, usually more setup is needed (see
    example for some details).
  }
  \item{fun}{
    if \code{fun} is NULL, the out-of-sample risk is returned.
    \code{fun}, as a function of \code{object}, may extract any other
    characteristic of the cross-validated models. These are returned as
    is.
  }
  \item{...}{
    additional arguments passed to \code{\link[parallel]{mclapply}}
    eventually.
  }
  \item{max}{
    vector of maximum values for each of the parameters. Per default
    \code{cvrisk} uses the current \code{mstop} values but any value
    might be specified.
  }
  \item{length.out}{
    number of grid points per parameter. If this is a scalar the same
    length is used for each parameter. For different grids per parameter
    use a vector. Per default 10 grid points are used per parameter.
  }
  \item{min}{
    minimum value of the grid (specify a vector for different values
    dependent on the parameter). Per default each grid starts at 1.
  }
  \item{log}{
    (logical) Should the grid be equidistant on the log scale (default)
    or on the linear scale?
  }
}
\details{
  The number of boosting iterations is a hyper-parameter of the boosting
  algorithms implemented in this package. Honest, i.e., cross-validated,
  estimates of the empirical risk for different stopping parameters
  \code{mstop} are computed by this function which can be utilized to
  choose an appropriate number of boosting iterations to be applied.
  For details see the definition of the generic function
  \code{\link[mboost]{cvrisk}} in package \pkg{mboost}.
}

\seealso{
  \code{\link{glmboostLSS}}, \code{\link{gamboostLSS}} and \code{\link{blackboostLSS}} for fitting
  of GAMLSS models.
}
 \examples{

## to be added

}
\keyword{methods}
