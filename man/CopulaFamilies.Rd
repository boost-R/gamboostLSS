\name{CopulaFamilies}
\alias{CopulaFamilies}
\alias{copula_tau}

\title{ Bivariate Copula Families for Boosting GAMLSS-models}

\description{
  Constructs a \code{\link{Families}} object for a bivariate response whose
  dependence structure is modelled by a copula, with each marginal
  distribution and the copula parameter receiving its own additive predictor.
}

\usage{
CopulaFamilies(marginal1, marginal2, copula = "gaussian", theta = 1,
               stabilization = c("none", "MAD", "L2"))
}
copula_tau(model, newdata = NULL)

\arguments{
  \item{marginal1, marginal2}{
    \code{\link{Families}} objects specifying the marginal distributions of
    the first and second response, e.g. \code{\link{GaussianLSS}()}. See
    Details for the supported combinations.
  }
  \item{copula}{
    character, one of \code{"gaussian"}, \code{"clayton"}, \code{"gumbel"}
    or \code{"frank"}.
  }
  \item{theta}{
    initial value of the copula dependence parameter, on the copula's own
    scale. Note this is only used until the offset is computed: at the start
    of the fit the offset for \code{theta} is derived from the empirical
    Kendall's tau of the responses where that falls inside the copula's
    representable range (see Details), which overrides this argument.
  }
  \item{stabilization}{
    governs if the negative gradient should be standardized in each boosting
    step. It can be either \code{"none"}, \code{"MAD"} or \code{"L2"}. See
    \code{\link{Families}}.
  }
  \item{model}{
  a fitted \code{mboostLSS} model whose families were built by 
  \code{CopulaFamilies}.
  }
  \item{newdata}{
  optionally, a data frame at which to evaluate the copula parameter. If omitted,
  the traning data is used.
  }
}

\details{
  For two continuous responses, Sklar's theorem gives the joint density as
  the product of the copula density evaluated at the marginal CDFs and the
  two marginal densities, so the joint log-likelihood decomposes into
  \deqn{\log c(F_1(y_1), F_2(y_2); \theta) + \log f_1(y_1) + \log f_2(y_2).}
  When a margin is discrete or binary, no such additive decomposition
  exists and the likelihood is instead built from differences of the copula
  CDF (discrete-discrete) or of its h-function (mixed cases).

  The combination of marginal types is detected automatically. Supported
  types are \emph{continuous} (Gaussian, Gamma, log-normal, Beta,
  Student t, log-logistic, Weibull), \emph{discrete} (negative binomial,
  ZIP, ZINB) and \emph{binary} (\code{\link{BernoulliLSS}}, logit, probit or 
  cloglog link).

  Each copula parameter is mapped to the real line by its own link:
  \code{tanh} for the Gaussian copula (\eqn{\rho \in (-1,1)}), \code{exp}
  for Clayton (\eqn{\theta > 0}), \code{exp(f) + 1} for Gumbel
  (\eqn{\theta \ge 1}) and the identity for Frank (\eqn{\theta \ne 0}).
  Predictions for \code{theta} obtained with \code{type = "response"} are
  therefore already on the copula's scale.

  The starting value for \code{theta} is derived from the empirical
  Kendall's tau of the two responses, inverted through the copula's
  tau-to-theta relationship. If the empirical tau lies outside what the
  copula can represent -- Clayton and Gumbel cannot express negative
  dependence -- a copula-specific static default is used instead.
  
  \code{copula_tau} converts the fitted copula parameter into Kendall's tau, 
  which unlike \code{theta} is comparable across copula families: a Clayton 
  \code{theta} of 2 and a Gumbel \code{theta} of 2 both correspond to a tau of 
  0.5. Not that Clayton and Gumbel cannot represent negative dependence, so a 
  positive tau from those families fitted to negatively associated data 
  indicates misspecification rather than a real positive association.
}

\value{
  An object of class \code{families}, with one component per distribution
  parameter: the parameters of \code{marginal1} suffixed with \code{1},
  those of \code{marginal2} suffixed with \code{2}, and \code{theta} for the
  copula. For two \code{\link{GaussianLSS}} margins the components are
  \code{mu1}, \code{sigma1}, \code{mu2}, \code{sigma2} and \code{theta}.

  The result can be passed to \code{\link{gamboostLSS}},
  \code{\link{glmboostLSS}} or \code{\link{gamboostLSS}} exactly like any
  other \code{families} object. The response must be a two-column matrix.
}

\references{
}

\seealso{
  \code{\link{Families}} for the marginal distributions,
  \code{\link{BernoulliLSS}} for binary margins,
  \code{\link{mboostLSS}} for fitting, and
  \code{\link{make_copula_grid}} for selecting \code{mstop} via
  \code{\link{cvrisk}}.
}

\examples{
set.seed(1)
n <- 100
x <- rnorm(n)
eps <- matrix(rnorm(2 * n), n, 2) \%*\% chol(matrix(c(1, .5, .5, 1), 2, 2))
y <- cbind(x + eps[, 1], -x + eps[, 2])

## Gaussian copula with two Gaussian margins
cf <- CopulaFamilies(GaussianLSS(), GaussianLSS(), copula = "gaussian")
names(cf)

fit <- gamboostLSS(y ~ bbs(x), families = cf, data = data.frame(x = x),
                   control = boost_control(mstop = 50, nu = 0.1))
                   
## type = "response" puts every parameter on its natural scale;
## for theta that is the copula's own scale (here: Pearson's rho)
pr <- predict(fit, newdata = data.frame(x = c(-1, 0, 1)), type = "response")
pr$theta

## a mixed binary-continuous response
y_mixed <- cbind(rbinom(n, 1, plogis(x)), rnorm(n))
cf_mixed <- CopulaFamilies(BernoulliLSS(), GaussianLSS(), copula = "clayton")
names(cf_mixed)
}

\keyword{ models }
\keyword{ distributions }