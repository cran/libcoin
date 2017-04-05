
\name{NEWS}
\title{NEWS file for the \pkg{libcoin} package}

\section{Changes in Version 0.9-2 (2017-04-04)}{
  \itemize{
    \item Fix protect problem in \code{R_MaximumTest} reported by Tomas Kalibera.
  }
}

\section{Changes in Version 0.9-1 (2017-02-06)}{
  \itemize{
    \item Calling \pkg{libcoin}'s C routines (using \code{.Call()}) now requires
          the entry points to be specified as \R objects, i.e., the use of character
          strings is no longer allowed.
    \item Stop when all observations are missing.
    \item 2d case with missings in at least one variable used incorrect
          sample sizes and thus produced wrong results.
  }
}

\section{Changes in Version 0.9-0 (2016-12-09)}{
  \itemize{
    \item \pkg{libcoin} published on CRAN.
  }
}

