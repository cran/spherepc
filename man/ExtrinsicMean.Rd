\name{ExtrinsicMean}
\alias{ExtrinsicMean}
\title{Finding Extrinsic Mean}
\description{
This function identifies the extrinsic mean of data on the unit 2-sphere.
}
\usage{
ExtrinsicMean(data, weights = rep(1, nrow(data)))
}
\arguments{
  \item{data}{matrix or data frame consisting of spatial locations with two columns.
  Each row represents longitude and latitude (denoted by degrees).}
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
it is possible that the extrinsic mean of the data set does not exists; for example,
a great circle.
}
\author{
Jongmin Lee
}
\references{
Jongmin Lee, Jang-Hyun Kim and Hee-Seok Oh. (2021). Spherical Principal Curves. IEEE Transactions on Pattern Analysis and Machine Intelligence, 43, 2165-2171. <doi.org/10.1109/TPAMI.2020.3025327>.
Jang-Hyun Kim, Jongmin Lee and Hee-Seok Oh. (2020). Spherical Principal Curves <arXiv:2003.02578>.
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
sphereplot::rgl.sphgrid()    
sphereplot::rgl.sphpoints(data, radius = 1, col = "blue", size = 12)
sphereplot::rgl.sphpoints(In.mean[1], In.mean[2], radius = 1, col = "red", size = 12)
sphereplot::rgl.sphpoints(Ex.mean[1], Ex.mean[2], radius = 1, col = "green", size = 12)
}
