
\name{LinStatExpCov}
\alias{LinStatExpCov}
\alias{lmult}
\title{
  Linear Statistics with Expectation and Covariance
}
\description{
  Strasser-Weber type linear statistics and their expectation
  and covariance under the independence hypothesis
}
\usage{
LinStatExpCov(X, Y, ix = NULL, iy = NULL, weights = integer(0),
              subset = integer(0), block = integer(0), checkNAs = TRUE,
              varonly = FALSE, nresample = 0, standardise = FALSE,
              tol = sqrt(.Machine$double.eps))
lmult(x, object)
}
\arguments{
  \item{X}{numeric matrix of transformations.}
  \item{Y}{numeric matrix of influence functions.}
  \item{ix}{an optional integer vector expanding \code{X}.}
  \item{iy}{an optional integer vector expanding \code{Y}.}
  \item{weights}{an optional integer vector of non-negative case weights.}
  \item{subset}{an optional integer vector defining a subset of observations.}
  \item{block}{an optional factor defining independent blocks of observations.}
  \item{checkNAs}{a logical for switching off missing value checks.  This
    included switching off checks for suitable values of \code{subset}.
    Use at your own risk.}
  \item{varonly}{a logical asking for variances only.}
  \item{nresample}{an integer defining the number of permuted statistics to draw.}
  \item{standardise}{a logical asking to standardise the permuted statistics.}
  \item{tol}{tolerance for zero variances.}
  \item{x}{a contrast matrix to be left-multiplied in case \code{X} was a factor.}
  \item{object}{an object of class \code{LinStatExpCov}.}
}
\details{
  The function, after minimal preprocessing, calls the underlying C code
  and computes the linear statistic, its expectation and covariance and,
  optionally, \code{nresample} samples from its permutation distribution.

  When both \code{ix} and \code{iy} are missing, the number of rows of
  \code{X} and \code{Y} is the same, ie the number of observations.

  When \code{X} is missing and \code{ix} a factor, the code proceeds as
  if \code{X} were a dummy matrix of \code{ix} without explicitly
  computing this matrix.

  Both \code{ix} and \code{iy} being present means the code treats them
  as subsetting vectors for \code{X} and \code{Y}.  Note that \code{ix = 0}
  or \code{iy = 0} means that the corresponding observation is missing
  and the first row or \code{X} and \code{Y} must be zero.

  \code{lmult} allows left-multiplication of a contrast matrix when \code{X}
  was (equivalent to) a factor.
}
\value{
  A list.
}
\references{
  Strasser, H. and Weber, C.  (1999).  On the asymptotic theory of permutation
  statistics.  \emph{Mathematical Methods of Statistics} \bold{8}(2), 220--250.
}
\examples{
wilcox.test(Ozone ~ Month, data = airquality, subset = Month \%in\% c(5, 8))

aq <- subset(airquality, Month \%in\% c(5, 8))
X <- as.double(aq$Month == 5)
Y <- as.double(rank(aq$Ozone))
doTest(LinStatExpCov(X, Y))
}
\keyword{htest}
