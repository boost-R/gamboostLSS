\name{cvrisk.mboostLSS}
\alias{cvrisk}
\alias{cvrisk.mboostLSS}
\alias{make.grid}

\title{ Cross-Validation }
\description{
  Multidimensional cross-validated estimation of the empirical risk for
  hyper-parameter selection.
}
\usage{
\method{cvrisk}{mboostLSS}(object, folds = cv(model.weights(object)),
       grid = make.grid(mstop(object)), papply = mclapply, fun = NULL, ...)
make.grid(max, length.out = 10, min = NULL, log = TRUE)
}

\arguments{
  \item{object}{
    an object of class \code{mboostLSS}, i.e. a boosted GAMLSS model.
  }
  \item{folds}{
    a weight matrix with number of rows equal to the number of
    observations. The number of columns corresponds to the number of
    cross-validation runs. Can be computed using function
    \code{\link{cv}} from package \pkg{mboost} and defaults to 25
    bootstrap samples.
  }
  \item{grid}{
    a list of vectors of stopping parameters the empirical risk is to be
    evaluated for. The number of vectors in the list must be equal to
    the number of parameters of the GAMLSS family. Per default,
    \code{make.grid(mstop(object))} is used.
  }
  \item{papply}{
    (parallel) apply function, defaults to  \code{\link[parallel]{mclapply}}.
    Alternatively, \code{\link[parallel]{parLapply}} can be used. In the
    latter case, usually more setup is needed (see example for some
    details). To run \code{cvrisk} sequentially (i.e. not in parallel),
    one can use \code{\link{lapply}}.
  }
  \item{fun}{
    if \code{fun} is NULL, the out-of-sample risk is returned. \code{fun},
    as a function of \code{object}, may extract any other characteristic
    of the cross-validated models. These are returned as is.
  }
  \item{\dots}{
    additional arguments passed to \code{\link[parallel]{mclapply}}.
  }
  \item{max}{
    a named vector of length equal to the number of parameters of the GAMLSS
    family (and names equal to the names of \code{families}) that
    determines the maximal values of the grid.
  }
  \item{length.out}{
    the number of grid points (default: 10). This can be either a vector
    of the same length as \code{max} (with different values) or a scalar
    (which is then used as length for all grids).
  }
  \item{min}{
    minimal value of the grid. Per default the grid starts at 1 but
    other values (smaller \code{max}) are possible. This can be either a
    vector of the same length as \code{max} (with different values) or a
    scalar (which is then used as \code{min} for all grids).
  }
  \item{log}{
    should the grid be on a logarithmic scale? Default: TRUE.
  }
}
\details{
  The number of boosting iterations is a hyper-parameter of the
  boosting algorithms implemented in this package. Honest,
  i.e., cross-validated, estimates of the empirical risk
  for different stopping parameters \code{mstop} are computed by
  this function which can be utilized to choose an appropriate
  number of boosting iterations to be applied. For details see
  \code{\link{cvrisk.mboost}}.

  \code{make.grid} eases the creation of an equidistand, integer-valued
  grids, which can be used with \code{cvrisk}. Per default, the grid is
  equidistant on a logarithmic scale.
}
\value{
  An object of class \code{cvriskLSS} (when \code{fun} wasn't
  specified), basically a matrix containing estimates of the empirical
  risk for a varying number of bootstrap iterations. \code{plot} and
  \code{print} methods are available as well as an \code{mstop} method.
}
\seealso{
  \code{\link{cvrisk.mboost}} and \code{\link{cv}} (both in package
  \pkg{mboost})
}
\examples{
## Data generating process:
set.seed(1907)
x1 <- rnorm(1000)
x2 <- rnorm(1000)
x3 <- rnorm(1000)
x4 <- rnorm(1000)
x5 <- rnorm(1000)
x6 <- rnorm(1000)
mu    <- exp(1.5 +1 * x1 +0.5 * x2 -0.5 * x3 -1 * x4)
sigma <- exp(-0.4 * x3 -0.2 * x4 +0.2 * x5 +0.4 * x6)
y <- numeric(1000)
for( i in 1:1000)
    y[i] <- rnbinom(1, size = sigma[i], mu = mu[i])
dat <- data.frame(x1, x2, x3, x4, x5, x6, y)

## linear model with y ~ . for both components: 400 boosting iterations
model <- glmboostLSS(y ~ ., families = NBinomialLSS(), data = dat,
                     control = boost_control(mstop = 100),
                     center = TRUE)

## set up a grid
grid <-  make.grid(mstop(model), length.out = 5)
grid
## a tiny toy example (5-fold bootsrap with maximum stopping value 100)
cvr <- cvrisk(model, folds = cv(model.weights(model), B = 5),
              grid = grid)
cvr
plot(cvr)
mstop(cvr)

## a more realistic example
\dontrun{  ## (automatically) as this really (!) takes a while
  ## set up a grid with maximum 400 for both components
  grid <- make.grid(c(mu = 400, sigma = 400))
  grid
  cvr <- cvrisk(model, grid = grid)
  mstop(cvr)
  ## reset model to mstop values:
  model[mstop(cvr)]
}
## Other grids:
make.grid(mstop(model), length.out = 3)
make.grid(c(mu = 400, sigma = 400), log = FALSE)
make.grid(c(mu = 400, sigma = 400), length.out = 4, min = 100, log = FALSE)
}
\keyword{models}
\keyword{regression}
