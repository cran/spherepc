\name{Trans.sph}
\alias{Trans.sph}
\title{Transforming into spherical coordinate}
\description{
This function converts a Euclidean coordinate to a spherical coordinate.
}
\usage{
Trans.sph(vec)
}
\arguments{
  \item{vec}{three-dimensional Euclidean coordinate.}
}
\details{
This function converts a three-dimensional Euclidean coordinate to a two-dimensional spherical coordinate.
If \code{vec} is not in the unit sphere, it is divided by its magnitude so that the result lies on the unit sphere.
}
\value{
two-dimensional vector.
}
\author{
Jongmin Lee
}
\seealso{
\code{\link{Trans.Euclid}}.
}
\examples{
Trans.sph(c(1, 0, 0))
Trans.sph(c(0, 1, 0))
Trans.sph(c(0, 0, 1))
Trans.sph(c(-1, 0 , 0))
Trans.sph(c(0, -1, 0))
}
