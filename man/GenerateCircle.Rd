\name{GenerateCircle}
\alias{GenerateCircle}
\title{Generating circle on sphere}
\description{
This function makes a circle on the sphere.
}
\usage{
GenerateCircle(center, radius, T = 1000)
}
\arguments{
  \item{center}{center of circle.}
  \item{radius}{radius of circle. It should be in [0, pi].}
  \item{T}{the number of points in circle.}
}
\details{
This function makes a circle on a sphere.
}
\value{
matrix consisting of spatial locations with two columns.
}
\author{
Jongmin Lee
}
\seealso{
\code{\link{PrincipalCircle}}.
}
\examples{
library(rgl)
library(sphereplot)
library(geosphere)
circle <- GenerateCircle(c(0, 0), 1)
# plot (It requires to load 'rgl', 'sphereplot', and 'geosphere' R package.)
rgl.sphgrid()
rgl.sphpoints(circle[, 1], circle[, 2], radius = 1, col = "blue", size = 12)
}
