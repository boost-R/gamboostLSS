gamboostLSS
===========

[![Build Status (Linux)](https://travis-ci.org/hofnerb/gamboostLSS.svg?branch=master)](https://travis-ci.org/hofnerb/gamboostLSS) 
[![Build status (Windows)](https://ci.appveyor.com/api/projects/status/81eo6c6v7v4h2llo/branch/master?svg=true)](https://ci.appveyor.com/project/hofnerb/gamboostlss/branch/master)
[![CRAN Status Badge](http://www.r-pkg.org/badges/version/gamboostLSS)](http://cran.r-project.org/package=gamboostLSS)
[![Coverage Status](https://coveralls.io/repos/hofnerb/gamboostLSS/badge.svg?branch=master&service=github)](https://coveralls.io/github/hofnerb/gamboostLSS?branch=master)
[![](http://cranlogs.r-pkg.org/badges/gamboostLSS)](http://cran.rstudio.com/web/packages/gamboostLSS/index.html)

`gamboostLSS` implements boosting algorithms for fitting generalized linear,
additive and interaction models for to potentially high-dimensional data.
Instead of modeling only the mean, `gamboostLSS` enables the user to model
various distribution parameters such as location, scale and shape at the same
time (hence the name GAMLSS, generalized additive models for location, scale and
shape).

## Installation

- Current version (from CRAN): 
  ```
  install.packages("gamboostLSS")
  ```

- Latest **patch version** (under development) from GitHub:
  ```
  library("devtools")
  install_github("boost-R/gamboostLSS/patch")
  library("gamboostLSS")
  ```

- Latest **development version** from GitHub:
  ```
  library("devtools")
  install_github("boost-R/gamboostLSS/pkg")
  library("gamboostLSS")
  ```

  To be able to use the `install_github()` command, one needs to install `devtools` first:
  ```
  install.packages("devtools")
  ```

## Using gamboostLSS

Instructions on how to use `gamboostLSS` can be found here:
- [gamboostLSS tutorial](http://arxiv.org/pdf/1407.1774v1); This is apreliminary version currently under review.
