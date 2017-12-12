
\name{ctabs}
\alias{ctabs}
\title{
    Cross Tabulation
}
\description{
    Efficient weighted cross tabulation of two factors and a block
}
\usage{
ctabs(ix, iy = integer(0), block = integer(0), weights = integer(0),
        subset = integer(0), checkNAs = TRUE)
}
\arguments{
  \item{ix}{a integer of positive values with zero indicating a missing.}
  \item{iy}{an optional integer of positive values with zero indicating a missing.}
  \item{block}{an optional blocking factor without missings.}
  \item{weights}{an optional vector of weights, integer or double.}
  \item{subset}{an optional integer vector indicating a subset.}
  \item{checkNAs}{a logical for switching off missing value checks.}
}
\details{
  A faster version of \code{xtabs(weights ~ ix + iy + block, subset)}.
}
\value{
  If \code{block} is present, a three-way table. Otherwise,
  a one- or two-dimensional table.
}
\examples{

  ctabs(ix = 1:5, iy = 1:5, weights = 1:5 / 5)

}
\keyword{univar}
