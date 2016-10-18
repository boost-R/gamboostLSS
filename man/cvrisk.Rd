\name{cvrisk.mboostLSS}
\alias{cvrisk}
\alias{cvrisk.mboostLSS}
\alias{make.grid}
\alias{plot.cvriskLSS}

\title{ Cross-Validation }
\description{
  Multidimensional cross-validated estimation of the empirical risk for
  hyper-parameter selection.
}
\usage{
\method{cvrisk}{mboostLSS}(object, folds = cv(model.weights(object)),
       grid = make.grid(mstop(object)), papply = mclapply,
       trace = TRUE, fun = NULL, ...)

make.grid(max, length.out = 10, min = NULL, log = TRUE,
          dense_mu_grid = TRUE)

\method{plot}{cvriskLSS}(x, type = c("heatmap", "lines"),
     xlab = NULL, ylab = NULL, ylim = range(x),
     main = attr(x, "type"), ...)
}

\arguments{
  \item{object}{
    an object of class \code{mboostLSS}, i.e., a boosted GAMLSS model.
  }
  \item{folds}{
    a weight matrix with number of rows equal to the number of
    observations. The number of columns corresponds to the number of
    cross-validation runs. Can be computed using function
    \code{\link{cv}} from package \pkg{mboost} and defaults to 25
    bootstrap samples.
  }
  \item{grid}{
    a matrix of stopping parameters the empirical risk is to be
    evaluated for. Each row represents a parameter combination. The
    number of columns must be equal to the number of parameters of the
    GAMLSS family. Per default, \code{make.grid(mstop(object))} is used.
  }
  \item{papply}{
    (parallel) apply function, defaults to  \code{\link[parallel]{mclapply}}.
    Alternatively, \code{\link[parallel]{parLapply}} can be used. In the
    latter case, usually more setup is needed. To run \code{cvrisk}
    sequentially (i.e. not in parallel), one can use \code{\link{lapply}}.
  }
  \item{trace}{
    should status information beein printed during cross-validation?
    Default: \code{TRUE}.
  }
  \item{fun}{
    if \code{fun} is NULL, the out-of-sample risk is returned. \code{fun},
    as a function of \code{object}, may extract any other characteristic
    of the cross-validated models. These are returned as is.
  }
  \item{\dots}{
    additional arguments passed to \code{\link[parallel]{mclapply}} or
    the \code{plot} function depending on the context.
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
    should the grid be on a logarithmic scale? Default: \code{TRUE}.
  }
  \item{dense_mu_grid}{
    should the grid in the \code{mu} component be extended for all
    values of the \code{mstop} values corresponding to \code{mu} that
    are greater or equal to all other parameters in this combination.
    These values can be computed without or with very little additional
    computational costs. For details see examples.
  }
  \item{x}{
    an object of class \code{cvriskLSS}, which results from running
    \code{cvrisk}.
  }
  \item{type}{
    should \code{"lines"} or a \code{"heatmap"} (default) be plotted?
    See details.
  }
  \item{xlab, ylab}{
    user-specified labels for the x-axis and y-axis of the plot (which
    are usually not needed). The defaults depend on the plot \code{type}.
  }
  \item{ylim}{
    limits of the y-axis. Only applicable for the line plot.
  }
  \item{main}{
    a title for the plots.
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

  The line plot depicts the avarage risk for each grid point and
  additionally shows information on the variability of the risk from
  fold to fold. The heatmap shows only the average risk but in a nicer
  fashion.

  Hofner et al. (2016) provide a detailed description of
  cross-validation for \code{\link{gamboostLSS}} models and show a
  worked example.
}
\value{
  An object of class \code{cvriskLSS} (when \code{fun} wasn't
  specified), basically a matrix containing estimates of the empirical
  risk for a varying number of bootstrap iterations. \code{plot} and
  \code{print} methods are available as well as an \code{mstop} method.
}
\references{
  B. Hofner, A. Mayr, M. Schmid (2016). gamboostLSS: An R Package for
  Model Building and Variable Selection in the GAMLSS Framework.
  Journal of Statistical Software, 74(1), 1-31.

  Available as \code{vignette("gamboostLSS_Tutorial")}.
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

## linear model with y ~ . for both components: 100 boosting iterations
model <- glmboostLSS(y ~ ., families = NBinomialLSS(), data = dat,
                     control = boost_control(mstop = 100),
                     center = TRUE)

## set up a grid
grid <-  make.grid(mstop(model), length.out = 5, dense_mu_grid = FALSE)
plot(grid)

\donttest{### Do not test the following code per default on CRAN as it takes some time to run:
### a tiny toy example (5-fold bootsrap with maximum stopping value 100)
## (to run it on multiple cores of a Linux or Mac OS computer remove
##  set papply = mclapply (default) and set mc.nodes to the
##  appropriate number of nodes)
cvr <- cvrisk(model, folds = cv(model.weights(model), B = 5),
              papply = lapply, grid = grid)
cvr
## plot the results
par(mfrow = c(1, 2))
plot(cvr)
plot(cvr, type = "lines")
## extract optimal mstop (here: grid to small)
mstop(cvr)
### END (don't test automatically)
}

\donttest{### Do not test the following code per default on CRAN as it takes some time to run:
### a more realistic example
grid <- make.grid(c(mu = 400, sigma = 400), dense_mu_grid = FALSE)
plot(grid)
cvr <- cvrisk(model, grid = grid)
mstop(cvr)
## set model to optimal values:
mstop(model) <- mstop(cvr)
### END (don't test automatically)
}

### Other grids:
plot(make.grid(mstop(model), length.out = 3, dense_mu_grid = FALSE))
plot(make.grid(c(mu = 400, sigma = 400), log = FALSE, dense_mu_grid = FALSE))
plot(make.grid(c(mu = 400, sigma = 400), length.out = 4,
               min = 100, log = FALSE, dense_mu_grid = FALSE))


### Now use dense mu grids
# standard grid
plot(make.grid(c(mu = 100, sigma = 100), dense = FALSE),
     pch = 20, col = "red")
# dense grid for all mstop_mu values greater than mstop_sigma
grid <- make.grid(c(mu = 100, sigma = 100))
points(grid, pch = 20, cex = 0.2)
abline(0,1)

# now with three parameters
grid <- make.grid(c(mu = 100, sigma = 100, df = 30),
                  length.out = c(5, 5, 2), dense = FALSE)
densegrid <- make.grid(c(mu = 100, sigma = 100, df = 30),
                       length.out = c(5, 5, 2))
par(mfrow = c(1,2))
# first for df = 1
plot(grid[grid$df == 1, 1:2], main = "df = 1", pch = 20, col = "red")
abline(0,1)
abline(v = 1)
# now expand grid for all mu values greater the corresponding sigma
# value (i.e. below the bisecting line) and above df (i.e. 1)
points(densegrid[densegrid$df == 1, 1:2], pch = 20, cex = 0.2)

# now for df = 30
plot(grid[grid$df == 30, 1:2], main = "df = 30", pch = 20, col = "red")
abline(0,1)
abline(v = 30)
# now expand grid for all mu values greater the corresponding sigma
# value (i.e. below the bisecting line) and above df (i.e. 30)
points(densegrid[densegrid$df == 30, 1:2], pch = 20, cex = 0.2)
}
\keyword{models}
\keyword{regression}
