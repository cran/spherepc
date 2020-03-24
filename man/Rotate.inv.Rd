\name{Rotate.inv}
\alias{Rotate.inv}
\title{Rotating point on sphere}
\description{
Rotate a point on the unit sphere.
}
\usage{
 Rotate.inv(pt1, pt2)
}
\arguments{
  \item{pt1}{a spatial location.}
  \item{pt2}{a spatial location.}
}
\details{
This function rotates \code{pt2} to the extent that the spherical coordinate (0, 90) is rotated to \code{pt1}.
The function is the inverse of the \code{Rotate} function, and returns a point as a form of three-dimensional Euclidean coordinate.
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
\code{\link{Rotate}}.
}
\examples{
## If "pt1" is north pole (= (0, 90)), then \code{Rotate.inv} returns Euclidean coordinate of "pt2".
# It returns Euclidean coornate of spatial location (-100, 80).
Rotate.inv(c(0, 90), c(-100, 80)) 
# It converts spatial coordinate (-100, 80) to Euclidean coordinate.
Trans.Euclid(c(-100, 80))                                                   
}
