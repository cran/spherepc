\name{Trans.Euclid}
\alias{Trans.Euclid}
\title{Transforming into Euclidean coordinate}
\description{
This function converts a spherical coordinate to a Euclidean coordinate.
}
\usage{
Trans.Euclid(vec)
}
\arguments{
  \item{vec}{two-dimensional spherical coordinate.}
}
\details{
This function converts a two-dimensional spherical coordinate to a three-dimensional Euclidean coordinate.
Longitude should be range from -180 to 180 and latitude from -90 to 90.
}
\value{
three-dimensional vector.
}
\author{
Jongmin Lee
}
\seealso{
\code{\link{Trans.sph}}.
}
\examples{
Trans.Euclid(c(0, 0))
Trans.Euclid(c(0, 90))
Trans.Euclid(c(90, 0))
Trans.Euclid(c(180, 0))
Trans.Euclid(c(-90, 0))
}
