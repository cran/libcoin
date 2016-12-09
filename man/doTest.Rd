\name{doTest}
\alias{doTest}
\title{
  Permutation Test
}
\description{
  Perform permutation test for a linear statistic
}
\usage{
doTest(object, teststat = c("maximum", "quadratic", "scalar"),
       alternative = c("two.sided", "less", "greater"),
       pvalue = TRUE, lower = FALSE, log = FALSE,
       minbucket = 10L, ordered = TRUE, pargs = GenzBretz())
}
\arguments{
  \item{object}{an object returned by \code{\link{LinStatExpCov}}.
}
  \item{teststat}{type of test statistic to use.
}
  \item{alternative}{alternative for scalar or maximum-type statistics.
}
  \item{pvalue}{a logical indicating if a p-value shall be computed.
}
  \item{lower}{a logical indicating if a p-value (\code{lower} is \code{FALSE})
               or 1 - p-value (\code{lower} is \code{TRUE}) shall be returned.
}
  \item{log}{a logical, if \code{TRUE} probabilities are log-probabilities.
}
  \item{minbucket}{minimum weight in either of two groups for maximally selected
                   statistics.
}
  \item{ordered}{a logical, if \code{TRUE} maximally selected statistics
                 assume that the cutpoints are ordered.
}
  \item{pargs}{arguments as in \code{\link[mvtnorm]{GenzBretz}}.
}
}
\details{
  Computes a test statistic, a corresponding p-value and, optionally, cutpoints for
  maximally selected statistics.
}
\value{
  A list.
}
\keyword{htest}
