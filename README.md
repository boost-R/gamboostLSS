gamboostLSS
===========

[![Build Status (Linux)](https://travis-ci.org/boost-R/gamboostLSS.svg?branch=devel)](https://travis-ci.org/boost-R/gamboostLSS)
[![Build status (Windows)](https://ci.appveyor.com/api/projects/status/373t0tvx5v1i5ooq/branch/devel?svg=true)](https://ci.appveyor.com/project/hofnerb/gamboostlss-s2whe/branch/devel)
[![CRAN Status Badge](http://www.r-pkg.org/badges/version/gamboostLSS)](http://cran.r-project.org/package=gamboostLSS)
[![Coverage Status](https://coveralls.io/repos/github/boost-R/gamboostLSS/badge.svg?branch=devel)](https://coveralls.io/github/boost-R/gamboostLSS?branch=devel)
[![](http://cranlogs.r-pkg.org/badges/gamboostLSS)](http://cran.rstudio.com/web/packages/gamboostLSS/index.html)


## Note

This is the version of `gamboostLSS` which was used for the Statistics & Computation submission.
It is currently not further developed, but required to keep our experiments reproducible.
Based on the results in the [paper]((https://arxiv.org/abs/1611.10171)), the `methode = "outer"` variant is not further developed.
Only use this version if you want to reproduce the experiments in the above mentioned paper.
Otherwise use the [developement](https://github.com/boost-R/gamboostLSS/tree/devel) or [stable](https://github.com/boost-R/gamboostLSS/tree/master) version of the package.



`gamboostLSS` implements boosting algorithms for fitting generalized linear,
additive and interaction models for to potentially high-dimensional data.
Instead of modeling only the mean, `gamboostLSS` enables the user to model
various distribution parameters such as location, scale and shape at the same
time (hence the name GAMLSS, generalized additive models for location, scale and
shape).


## Using gamboostLSS

For installation instructions see below.

Instructions on how to use `gamboostLSS` can be found here:
- [gamboostLSS tutorial](https://www.jstatsoft.org/article/view/v074i01).

Details on the noncyclical fitting method can be found here:
- [noncyclical fitting](https://arxiv.org/abs/1611.10171); This is a preleminary version currently under review.

## Issues & Feature Requests

For issues, bugs, feature requests etc. please use the [GitHub Issues](https://github.com/boost-R/gamboostLSS/issues).

## Installation

- Current version (from CRAN):
  ```
  install.packages("gamboostLSS")
  ```

- Latest **patch version** (patched version of CRAN package; under development) from GitHub:
  ```
  library("devtools")
  install_github("boost-R/gamboostLSS")
  library("gamboostLSS")
  ```

- Latest **development version** (version with new features; under development) from GitHub:
  ```
  library("devtools")
  install_github("boost-R/gamboostLSS", ref = "devel")
  library("gamboostLSS")
  ```

  To be able to use the `install_github()` command, one needs to install `devtools` first:
  ```
  install.packages("devtools")
  ```

