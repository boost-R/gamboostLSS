\name{stabsel}
\alias{stabsel.mboostLSS}
\alias{selected.stabsel_mboostLSS}

\title{
    Stability Selection
}
\description{
    Selection of influential variables or model components with error control.
}
\usage{
## a method to compute stability selection paths for fitted mboostLSS models
\method{stabsel}{mboostLSS}(x, cutoff, q, PFER, mstop = NULL,
        folds = subsample(model.weights(x), B = B),
        B = ifelse(sampling.type == "MB", 100, 50),
        assumption = c("unimodal", "r-concave", "none"),
        sampling.type = c("SS", "MB"),
        papply = mclapply, verbose = TRUE, FWER, eval = TRUE, ...)
## a method to get the selected parameters
\method{selected}{stabsel_mboostLSS}(object, parameter = NULL, ...)
}
\arguments{
  \item{x}{an fitted model of class \code{"mboostLSS"} or \code{"nc_mboostLSS"}.}
  \item{cutoff}{cutoff between 0.5 and 1. Preferably a value between 0.6
    and 0.9 should be used.}
  \item{q}{number of (unique) selected variables (or groups of variables
    depending on the model) that are selected on each subsample.}
  \item{PFER}{upper bound for the per-family error rate. This
    specifies the amount of falsely selected base-learners, which is
    tolerated. See details.}
  \item{mstop}{mstop value to use, if no value is supplied the mstop value of the fitted model is used.}
  \item{folds}{ a weight matrix with number of rows equal to the number
    of observations, see \code{\link{cvrisk}} and
    \code{\link[stabs]{subsample}}. Usually one should not
    change the default here as subsampling with a fraction of \eqn{1/2}
    is needed for the error bounds to hold. One usage scenario where
    specifying the folds by hand might be the case when one has
    dependent data (e.g. clusters) and thus wants to draw clusters
    (i.e., multiple rows together) not individuals.}
  \item{assumption}{ Defines the type of assumptions on the
    distributions of the selection probabilities and simultaneous
    selection probabilities. Only applicable for
    \code{sampling.type = "SS"}. For \code{sampling.type = "MB"} we
    always use code{"none"}.}
  \item{sampling.type}{ use sampling scheme of of Shah & Samworth
    (2013), i.e., with complementarty pairs (\code{sampling.type = "SS"}),
    or the original sampling scheme of Meinshausen & Buehlmann (2010).}
  \item{B}{ number of subsampling replicates. Per default, we use 50
    complementary pairs for the error bounds of Shah & Samworth (2013)
    and 100 for the error bound derived in  Meinshausen & Buehlmann
    (2010). As we use \eqn{B} complementray pairs in the former case
    this leads to \eqn{2B} subsamples.}
  \item{papply}{ (parallel) apply function, defaults to
    \code{\link[parallel]{mclapply}}. Alternatively, \code{parLapply}
    can be used. In the latter case, usually more setup is needed (see
    example of \code{\link{cvrisk}} for some details).}
  \item{verbose}{ logical (default: \code{TRUE}) that determines wether
    \code{warnings} should be issued. }
  \item{FWER}{ deprecated. Only for compatibility with older versions,
    use PFER instead.}
  \item{eval}{ logical. Determines whether stability selection is
    evaluated (\code{eval = TRUE}; default) or if only the parameter
    combination is returned.}
  \item{object}{ a object of class \code{"stabsel_mboostLSS"}.}
  \item{parameter}{ select one or multiple effects.}
  \item{\dots}{ additional arguments to parallel apply methods such as
    \code{\link{mclapply}} and to \code{\link{cvrisk}}.}
}
\details{

  For details see \code{\link[stabs]{stabsel}} in package \pkg{stabs}
  and Hofner et al. (2014).

}
\value{
  An object of class \code{stabsel} with a special \code{print} method.
  The object has the following elements:
  \item{phat}{selection probabilities.}
  \item{selected}{elements with maximal selection probability greater
    \code{cutoff}.}
  \item{max}{maximum of selection probabilities.}
  \item{cutoff}{cutoff used.}
  \item{q}{average number of selected variables used.}
  \item{PFER}{per-family error rate.}
  \item{sampling.type}{the sampling type used for stability selection.}
  \item{assumption}{the assumptions made on the selection
    probabilities.}
  \item{call}{the call.}
}
\references{

  B. Hofner, L. Boccuto and M. Goeker (2014),
  Controlling false discoveries in high-dimensional situations: Boosting
  with stability selection. \emph{Technical Report}, arXiv:1411.1285.\cr
  \url{http://arxiv.org/abs/1411.1285}.

  N. Meinshausen and P. Buehlmann (2010), Stability selection.
  \emph{Journal of the Royal Statistical Society, Series B},
  \bold{72}, 417--473.

  R.D. Shah and R.J. Samworth (2013), Variable selection with error
  control: another look at stability selection. \emph{Journal of the Royal
  Statistical Society, Series B}, \bold{75}, 55--80.

}
\seealso{
  \code{\link[stabs]{stabsel}} and
  \code{\link[stabs]{stabsel_parameters}}
}
\examples{

### Data generating process:
set.seed(1907)
x1 <- rnorm(500)
x2 <- rnorm(500)
x3 <- rnorm(500)
x4 <- rnorm(500)
x5 <- rnorm(500)
x6 <- rnorm(500)
mu    <- exp(1.5 +1 * x1 +0.5 * x2 -0.5 * x3 -1 * x4)
sigma <- exp(-0.4 * x3 -0.2 * x4 +0.2 * x5 +0.4 * x6)
y <- numeric(500)
for( i in 1:500)
    y[i] <- rnbinom(1, size = sigma[i], mu = mu[i])
dat <- data.frame(x1, x2, x3, x4, x5, x6, y)

### linear model with y ~ . for both components: 400 boosting iterations
model <- glmboostLSS(y ~ ., families = NBinomialLSS(), data = dat,
                     control = boost_control(mstop = 400),
                     center = TRUE, method = "noncyclic")

\donttest{### Do not test the following code per default on CRAN as it takes some time to run:

#run stability selection 
(s <- stabsel(model, q = 5, PFER = 1))
#get selected effects
selected(s)

#visualize selection frequencies 
plot(s)

### END (don't test automatically)
}
}
\keyword{nonparametric}