\name{ctabs}
\alias{ctabs}
\title{
    Cross Tabulation
}
\description{
    Efficient weighted cross tabulation of two factors and a block
}
\usage{
ctabs(ix, iy = integer(0), weights = integer(0),
      subset = integer(0), block = integer(0))
}
\arguments{
  \item{ix}{a integer of positive values with zero indicating a missing.
}
  \item{iy}{an optional integer of positive values with zero indicating a missing.
}
  \item{weights}{a vector of weights, integer or double.
}
  \item{subset}{an optional integer vector indicating a subset.
}
  \item{block}{an optional blocking factor without missings.
}
}
\details{
  A faster version of \code{xtabs(weights ~ ix + iy + block, subset)}.
}
\value{
  If \code{block} is present, a three-way table. Otherwise,
  a one- or two-dimensional table.
}
\examples{

  ctabs(1:5, 1:5, 1:5 / 5)

}
\keyword{univar}
