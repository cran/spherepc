\name{Rotate}
\alias{Rotate}
\title{Rotating point on a sphere}
\description{
Rotate a point on the unit 2-sphere.
}
\usage{
Rotate(pt1, pt2)
}
\arguments{
  \item{pt1}{spatial location.}
  \item{pt2}{spatial location.}
}
\details{
This function rotates \code{pt2} to the extent that \code{pt1} to spherical coordinate (0, 90).
The function returns a point as a form of three-dimensional Euclidean coordinate.
}
\value{
three-dimensional vector.
}
\references{
https://en.wikipedia.org/wiki/Rodrigues_rotation_formula
}
\author{
Jongmin Lee
}
\seealso{
\code{\link{Rotate.inv}}.
}
\examples{
## If "pt1" is north pole (= (0, 90)), Rotate() function returns Euclidean coordinate of "pt2".
Rotate(c(0, 90), c(10, 10))   # It returns Euclidean coordinate of spatial location (10, 10).
# The Trans.Euclid() function converts spatial coordinate (10, 10) to Euclidean coordinate.
Trans.Euclid(c(10, 10))       
}
