\name{PGA}
\alias{PGA}
\title{Principal geodesic analysis}
\description{
This function performs principal geodesic analysis.
}
\usage{
PGA(data, col = c("blue", "red"), size = c(12, 6))
}
\arguments{
  \item{data}{matrix or data frame consisting of spatial locations with two columns.
  Each row represents longitude and latitude.}
  \item{col}{two-dimensional vector which represents colors of data and the principal geodesic line.}
  \item{size}{two-dimensional vector which represents sizes of data and the principal geodesic line.}
}
\details{
This function performs principal geodesic analysis.
}
\value{
plot and a list consisting of
 \item{line}{spatial locations of points in the principal geodesic line.}
}
\references{
Fletcher, P. T., Lu, C., Pizer, S. M. and Joshi, S. (2004). Principal geodesic analysis for the study of nonlinear statistics of shape. IEEE Transactions on Medical Imaging, 23, 995-1005.
}
\author{
Jongmin Lee
}
\note{
This function requires to load 'sphereplot', 'geosphere' and 'rgl' R package.
}
\seealso{
\code{\link{LPG}}.
}
\examples{
library(rgl)
library(sphereplot)
library(geosphere)
#### example 1: noisy half-great circle data
circle <- GenerateCircle(c(150, 60), radius = pi / 2)
half.circle <- circle[circle[, 1] < 0, , drop = FALSE]
sigma <- 2
half.circle <- half.circle + sigma * rnorm(nrow(half.circle))
PGA(half.circle)

#### example 2: noisy S-shaped data 
#### The data consists of two parts: x ~ Uniform[0, 20], y = sqrt(20 * x - x^2) + N(0, sigma^2), 
#### x ~ Uniform[-20, 0], y = -sqrt(-20 * x - x^2) + N(0, sigma^2).
n <- 500                                           
x <- 60 * runif(n)
sigma <- 2                                         
y <- (60 * x - x^2)^(1/2) + sigma * rnorm(n)
simul.S1 <- cbind(x, y)
z <- -60 * runif(n)
w <- -(-60 * z - z^2)^(1/2)+ sigma * rnorm(n)
simul.S2 <- cbind(z, w)
simul.S <- rbind(simul.S1, simul.S2)
PGA(simul.S)
}
\keyword{ ~principal geodesic analysis }
