\name{india}
\alias{india}
\alias{india.bnd}
\docType{data}

\title{
  Malnutrition of Children in India (DHS, 1998-99)
}

\description{
  Data sample from the Standard Demographic and Health Survey, 1998-99,
  on malnutrition of children in India. The data set contains
  approximately 12\% of the observations in the original data set and
  only a (very small) subset  of variables. Additionally, a boundary
  file representing the districts of India is provided for spatial analysis.
}

\usage{
data(india)
data(india.bnd)
}
\details{
  For details on the boundary file see function \code{\link[BayesX]{read.bnd}}
  from package \pkg{BayesX}.
}

\format{
  A data frame with 4000 observations on the following 6 variables:
  \describe{
    \item{\code{stunting}}{A numeric score for malnutrition which ranges
      from -600 to 600, where negative values represent malnourished
      children.}
    \item{\code{cbmi}}{BMI of the child.}
    \item{\code{cage}}{Age of the child in months.}
    \item{\code{mbmi}}{BMI of the mother.}
    \item{\code{mage}}{Age of the mother in years.}
    \item{\code{mcdist}}{The district in India, where mother and child live.
      A factor encoded to match the map \code{\link{india.bnd}}.}
    \item{\code{mcdist_lab}}{The district in India, where mother and child live.
      A factor with actual district names.}
  }
}

\source{ The complete data set is provided by the Monitoring and
  Evaluation to Assess and Use Results Demographic and Health Surveys
  (MEASURE DHS) which is funded by the U.S. Agency of International
  Development (USAID). It can be obtained for research purposes (after
  registration) from
  \url{http://www.measuredhs.com/what-we-do/survey/survey-display-156.cfm}
  (Data set for All-India, Children's Recode: iakr42dt.zip) }
  \references{ For details on the data set see also:

  Fahrmeir L and Kneib T (2011),
  \emph{Bayesian smoothing and regression for longitudinal, spatial and
    event history data}, Oxford University Press.

}
\examples{
if (require("BayesX")) {
  ## plot distribution of stunting in India
  drawmap(india, map = india.bnd, regionvar = "mcdist", plotvar = "stunting")
}
}
\keyword{datasets}
