\name{ExtrinsicMean}
\alias{ExtrinsicMean}
\title{Finding Extrinsic Mean}
\description{
This function identifies the extrinsic mean of data on the sphere.
}
\usage{
ExtrinsicMean(data, weights = rep(1, nrow(data)))
}
\arguments{
  \item{data}{matrix or data frame consisting of spatial locations with two columns.
  Each row represents longitude and latitude.}
  \item{weights}{vector of weights.}
}
\details{
This function identifies the extrinsic mean of data.
}
\value{
two-dimensional vector.
}
\note{
In the case of spheres, if data set is not contained in a hemisphere, then 
it is possible that the extrinsic mean of the data set does not exists, for example,
a great circle.
}
\author{
Jongmin Lee
}
\references{
Jang-Hyun Kim, Jongmin Lee, Hee-Seok Oh. (2020). Spherical Principal Curves <arXiv:2003.02578>.
}
\seealso{
\code{\link{IntrinsicMean}}.
}
\examples{
#### comparison of Intrinsic mean and extrinsic mean.
#### example: noisy circular data set.
library(rgl)
library(sphereplot)
library(geosphere)
n <- 500                  # the number of samples.
x <- 360 * runif(n) - 180
sigma <- 5
y <- 60 + sigma * rnorm(n)
simul.circle <- cbind(x, y)
data <- simul.circle
In.mean <- IntrinsicMean(data)
Ex.mean <- ExtrinsicMean(data)
## plot (color of data is "blue"; that of intrinsic mean is "red" and 
## that of extrinsic mean is "green".)
rgl.sphgrid()    
rgl.sphpoints(data, radius = 1, col = "blue", size = 12)
rgl.sphpoints(In.mean[1], In.mean[2], radius = 1, col = "red", size = 12)
rgl.sphpoints(Ex.mean[1], Ex.mean[2], radius = 1, col = "green", size = 12)
}
