\name{IntrinsicMean}
\alias{IntrinsicMean}
\title{Finding Intrinsic Mean}
\description{
This function calculates the intrinsic mean of data on sphere.
}
\usage{
IntrinsicMean(data, weights = rep(1, nrow(data)), thres = 1e-5)
}
\arguments{
  \item{data}{matrix or data frame consisting of spatial locations with two columns.
  Each row represents longitude and latitude (denoted by degrees).}
  \item{weights}{vector of weights.}
  \item{thres}{threshold of the stopping conditions.}
}
\details{
This function calculates the intrinsic mean of data. The intrinsic mean is found by the gradient descent algorithm, which works well if the data is well-localized. In the case of spheres, if data is contained in a hemisphere, then the algorithm converges. 
}
\value{
two-dimensional vector.
}
\references{
Fletcher, P. T., Lu, C., Pizer, S. M. and Joshi, S. (2004). Principal geodesic analysis for the study of nonlinear statistics of shape. IEEE Transactions on Medical Imaging, 23, 995-1005.
Jongmin Lee, Jang-Hyun Kim and Hee-Seok Oh. (2021). Spherical Principal Curves. IEEE Transactions on Pattern Analysis and Machine Intelligence, 43. 2165-2171. <doi.org/10.1109/TPAMI.2020.3025327>.
}
\author{
Jongmin Lee
}
\seealso{
\code{\link{ExtrinsicMean}}.
}
\examples{
#### comparison of Intrinsic mean and extrinsic mean.
#### example: circular data set.
library(rgl)
library(sphereplot)
library(geosphere)
n <- 500
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
