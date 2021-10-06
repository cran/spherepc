\name{GenerateCircle}
\alias{GenerateCircle}
\title{Generating circle on sphere}
\description{
This function makes a circle on the unit 2-sphere.
}
\usage{
GenerateCircle(center, radius, T = 1000)
}
\arguments{
  \item{center}{center of circle with spatial locations (longitude and latitude denoted by degrees).}
  \item{radius}{radius of circle. It should be in [0, pi].}
  \item{T}{the number of points that make up circle. The points in circle are equally spaced. The default is 1000.}
}
\details{
This function makes a circle on the unit 2-sphere.
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
sphereplot::rgl.sphgrid()
sphereplot::rgl.sphpoints(circle[, 1], circle[, 2], radius = 1, col = "blue", size = 12)
}
\keyword{ ~exact principal circle }