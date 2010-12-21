\name{gamboostLSS-package}
\alias{gamboostLSS-package}

\docType{package}
\title{ 
  Boosting algorithms for GAMLSS
}
\description{
  Boosting methods for fitting generalized additive models for
  location, scale and shape (GAMLSS).
}
\details{

This package uses boosting algorithms for fitting GAMLSS (generalized
additive models for location, scale and shape). For information on
GAMLSS theory see Rigby and Stasinopoulos (2005), or the information
provided at \url{http://gamlss.org}.


The fitting methods \code{\link{glmboostLSS}} and
\code{\link{gamboostLSS}}, are alternatives for the algorithms provided
with \code{\link[gamlss]{gamlss}} in the \code{gamlss} package. They
offer shrinkage of effect estimates, intrinsic variable selecion and
model choice for potentially high-dimensional data settings.


\code{\link{glmboostLSS}} (for linear effects) and
\code{\link{gamboostLSS}} (for smooth effects) depend on their analogous
companions \code{\link[mboost]{glmboost}} and
\code{\link[mboost]{gamboost}} for generalized additive models
(contained in package \code{\link{mboost}}, see Hothorn et al.
2010a, b) and are similar in their usage.

The package includes some pre-defined GAMLSS distributions, but the user
can also specify new distributions with \code{\link{Families}}.


A wide range of different base-learners is available for covariate
effects (see \code{\link[mboost]{baselearners}}) including also random
(\code{brandom}) or spatial effects (\code{bspatial} or Markov random
fields \code{bmrf}). The selection of base-learnes is crucial as it
implies the kind of effect the covariate has on each distribution
parameter in the final GAMLSS.

}
\author{
Benjamin Hofner, Andreas Mayr, Nora Fenske, Matthias Schmid

Maintainer:  Benjamin Hofner <benjamin.hofner@imbe.med.uni-erlangen.de>
}
\references{

Mayr, A., Fenske, N., Hofner, B., Kneib, T. and Schmid, M. (2010): GAMLSS 
for high-dimensional data - a flexible approach based on boosting. 
Department of Statistics: Technical Report 98. \url{http://epub.ub.uni-muenchen.de/11938/}

Rigby, R. A. and D. M. Stasinopoulos (2005). Generalized additive models for location, scale and
shape (with discussion). Applied Statistics 54, 507-554.

Stasinopoulos, D. M. and R. A. Rigby (2007). Generalized additive models for location scale and
shape (GAMLSS) in R. Journal of Statistical Software 23(7).

Buehlmann, P. and Hothorn, T. (2007). Boosting algorithms: Regularization, 
prediction and model fitting. Statistical Science, 22(4), 477--505. 

Hothorn, T., Buehlmann, P., Kneib, T., Schmid, M. and Hofner, B. (2010a).
Model-based boosting 2.0. Journal of Machine Learning Research 11(Aug),
2109-2113.

Hothorn, T., Buehlmann, P., Kneib, T., Schmid, M. and Hofner, B. (2010b).
mboost: Model-based boosting 2.0. R package version 2.0-9.
\url{http://cran.r-project.org/web/packages/mboost}

}
\keyword{ package }
\seealso{ The \code{\link{mboost}} package for more on model-based boosting, or the 
          \code{\link[gamlss]{gamlss}} package for the original GAMLSS algorithms 
          provided by Rigby and Stasinopoulos.
}
\examples{

# Generate covariates
x1 <- runif(100)
x2 <- runif(100)
eta_mu <-     2 - 2*x1
eta_sigma <-  -1  + 2*x2

# Generate response: Negative Binomial Distribution
y <- numeric(100)
for( i in 1:100)  y[i] <- rnbinom(1, size=exp(eta_sigma[i]), mu=exp(eta_mu[i]))

# Model fitting, 300 boosting steps
mod1 <- glmboostLSS( y ~ x1 + x2, families=NBinomialLSS(), 
        control=boost_control(mstop=300))

# Shrinked effect estimates
coef(mod1, off2int=TRUE)

# Empirical risk with respect to mu
plot(risk(mod1)$mu)

# Empirical risk with respect to sigma
plot(risk(mod1)$sigma)
}
