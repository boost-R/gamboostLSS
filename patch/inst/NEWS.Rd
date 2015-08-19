\name{NEWS}
\title{News for Package 'gamboostLSS'}

\section{Changes in gamboostLSS version 1.2-0 (2015-08-19)}{
  \subsection{User-visible changes}{
    \itemize{
      \item New interface: Stabilization of negative gradient is now
      specified via families, e.g., \code{GaussianLSS(stabilization =
	"MAD")}. Using \code{options(gamboostLSS_stab_ngrad = TRUE)} is
      still supported for backward compatibility but discouraged.
      \item Added \code{weighted.median} which is used in \code{"MAD"}
      stabilization (see bugfixes).
      \item Added preliminary version of the tutorial paper as vignette,
      which is now available as \code{vignette("gamboostLSS_Tutorial",
	package = "gamboostLSS")}.
      \item Added \code{summary} function (request
      \href{https://github.com/hofnerb/gamboostLSS/issues/2}{#2}).
    }
  }
  \subsection{Miscellaneous}{
    \itemize{
      \item Added preliminary version of tutorial paper as reference.
    }
  }
  \subsection{Bug-fixes}{
    \itemize{
      \item Improve stored \code{call}s (closes
      \href{https://github.com/hofnerb/gamboostLSS/issues/3}{#3}).
      \item Use \code{weighted.median} to compute MAD when negative
      gradients are stabilized.
      \item Fixed handling of \pkg{gamlss.dist} families with \code{type
	= "Mixed"}.
      \item Fixed a bug in \code{as.families} with truncated \code{families}
      (which can be obtained using the package \pkg{gamlss.tr}).
      \item The internal helper \code{get_qfun} is now a generic
      function.
      \item  Adhere to CRAN policies regarding import of base packages
      (closes \href{https://github.com/hofnerb/gamboostLSS/issues/1}{#1}).
    }
  }
}

\section{Changes in gamboostLSS version 1.1-3 (2015-01-12)}{
  \subsection{Miscellaneous}{
    \itemize{
      \item Changed plain text \file{NEWS} to \file{inst/NEWS.Rd}
    }
  }
  \subsection{Bug-fixes}{
    \itemize{
      \item Changes in \file{inst/CITATION} to make CRAN happy:
      Citations can now be extracted without the need to install the
      package.
      \item Corrected \code{india} data set:\cr
      stunting ranged from -600 to 600 but should range from -6 to 6
      \item Fixed unstated dependencies
      \item Changed \code{require} to \code{requireNamespace}
      and added \code{gamlss.dist::} where necessary
      \item Fixed bug in \code{selected.mboostLSS}
      \item \code{predict()} returned a vector instead of a matrix
      \item \code{factor}s were coerced to \code{character}s in
      \code{predint()}
    }
  }
}


\section{Changes in gamboostLSS version 1.1-2 (2014-06-25)}{
  \subsection{Miscellaneous}{
    \itemize{
      \item added \pkg{survival} to suggests and load it in tests to
      reflect recent changes in \pkg{mboost}
      \item removed the generic function \code{risk()} as this function is
      now defined in \pkg{mboost}
    }
  }
  \subsection{Bug-fixes}{
    \itemize{
      \item check if response is correct for the \code{family} before
      the \code{offset} is computed
    }
  }
}

\section{Changes in gamboostLSS version 1.1-1 (2014-06-16)}{
  \subsection{Miscellaneous}{
    \itemize{
      \item re-added tests
      \item changed \code{\\dontrun} to \code{\\donttest} in examples
    }
  }
}

\section{Changes in gamboostLSS version 1.1-0 (2014-06-12)}{
  \subsection{User-visible changes}{
    \itemize{
       \item added new function \code{cvrisk} for generic
       cross-validation (\pkg{gamboostLSS} now requires \pkg{mboost} >=
       2.2-3 for the generic function definition)
       \item added a function to obtain cross-validation grids
       (\code{make_grid})
       \item added new function \code{as.families} to include
       \pkg{gamlss} families in \pkg{gamboostLSS}
       \item added new families \code{BetaLSS}, \code{GammaLSS},
       \code{GaussianLSS}, \code{ZIPoLSS}, \code{ZINBLSS}
       \item \code{GaussianLSS} is now the default for \code{families}
       \item added new functions \code{predint()} and
       \code{plot(predint())} to obtain and plot marginal prediction
       intervals
       \item added stabilization factor for negative gradients (MAD):
       turn on the stabilization by using\cr
       \code{options(gamboostLSS_stab_ngrad = TRUE)}
       \item added new data sets \code{india} and \code{india.bnd}
    }
  }
  \subsection{Miscellaneous}{
    \itemize{
      \item lots of minor improvements in the interface and methods
      \item improved manuals in various places
    }
  }
  \subsection{Bug-fixes}{
    \itemize{
      \item Speed of model fitting was dependent on iteration,
      i.e. the algorithm got slower for larger \code{mstop} values
    }
  }
}

\section{Changes in gamboostLSS version 1.0-3 (2011-11-16)}{
  \subsection{Miscellaneous}{
    \itemize{
      \item removed duplicated (and slightly altered) code from
      \pkg{mboost} that was required to make \pkg{gamboostLSS} work with
      earlier versions of \pkg{mboost};\cr
      \pkg{gamboostLSS} now requires \pkg{mboost} >= 2.0-12
    }
  }
  \subsection{Bug-fixes}{
    \itemize{
      \item fixed potential problem that might occur if users specify
      different outcomes for different components
      \item fixed bug in \file{tests/regtest-gamboostLSS.R}
    }
  }
}

\section{Changes in gamboostLSS version 1.0-2 (2011-08-31)}{
  \subsection{Bug-fixes}{
    \itemize{
      \item changed the way we use environments
      (to get rid of note: \code{no visible binding})
    }
  }
}

\section{Changes in gamboostLSS version 1.0-1 (2011-08-10)}{
  \subsection{Miscellaneous}{
    \itemize{
      \item updated dependencies
    }
  }
  \subsection{Bug-fixes}{
    \itemize{
      \item fixed bug in \code{selected()}
    }
  }
}

\section{Changes in gamboostLSS version 1.0 (2011-08-10)}{
  \subsection{Miscellaneous}{
    \itemize{
      \item updated references
      \item moved \pkg{gamboostLSS} to CRAN
    }
  }
  \subsection{Bug-fixes}{
    \itemize{
      \item some changes to make \pkg{gamboostLSS} work with
      \pkg{mboost} 2.0.x series
    }
  }
}

\section{Changes in gamboostLSS version 0.5 (2010-07-02)}{
  \itemize{
    \item initial development version of \pkg{gamboostLSS} on R-forge
  }
}
