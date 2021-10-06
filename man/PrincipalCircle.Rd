\name{PrincipalCircle}
\alias{PrincipalCircle}
\title{Principal circle on a sphere}
\description{
This function fits a principal circle on sphere via gradient descent algorithm.
}
\usage{
PrincipalCircle(data, step.size = 1e-3, thres = 1e-5, maxit = 10000)
}
\arguments{
  \item{data}{matrix or data frame consisting of spatial locations (longitude and latitude denoted by degrees) with two columns.}
  \item{step.size}{step size of gradient descent algorithm. For convergence of the algorithm, \code{step.size} is recommended to be below 0.01. The default is 1e-3.}
  \item{thres}{threshold of the stopping condition. The default is 1e-5.}
  \item{maxit}{maximum number of iterations. The default is 10000.}
}
\details{
This function fits a principal circle on sphere via gradient descent algorithm. The function returns three-dimensional vectors whose components represent longitude and latitude of the center and the radius of the circle in regular order.
}
\value{
three-dimensional vector.
}
\references{
Jongmin Lee, Jang-Hyun Kim and Hee-Seok Oh. (2021). Spherical principal curves. IEEE Transactions on Pattern Analysis and Machine Intelligence, 43, 2165-2171. <doi.org/10.1109/TPAMI.2020.3025327>.

Jongmin Lee, Jang-Hyun Kim and Hee-Seok Oh (2020), Spherical principal curves <arXiv:2003.02578>.
}
\author{
Jongmin Lee
}
\seealso{
\code{\link{GenerateCircle}}
}
\examples{
library(rgl)
library(sphereplot)
library(geosphere)
#### example 1: half-great circle data
circle <- GenerateCircle(c(150, 60), radius = pi/2)
half.great.circle <- circle[circle[, 1] < 0, , drop = FALSE]
sigma <- 2
half.great.circle <- half.great.circle + sigma * rnorm(nrow(half.great.circle))
## find a principal circle
\donttest{
PC <- PrincipalCircle(half.great.circle)
result <- GenerateCircle(PC[1:2], PC[3])
## plot
rgl.sphgrid()
rgl.sphpoints(half.great.circle, radius = 1, col = "blue", size = 12)
rgl.sphpoints(result, radius = 1, col = "red", size = 6)

#### example 2: circular data
n <- 700                            
x <- seq(-180, 180, length.out = n)
sigma <- 5                          
y <- 45 + sigma * rnorm(n)
simul.circle <- cbind(x, y)
## find a principal circle
PC <- PrincipalCircle(simul.circle)
result <- GenerateCircle(PC[1:2], PC[3])
## plot
sphereplot::rgl.sphgrid()
sphereplot::rgl.sphpoints(simul.circle, radius = 1, col = "blue", size = 12)
sphereplot::rgl.sphpoints(result, radius = 1, col = "red", size = 6)

#### example 3: earthquake data
data(Earthquake)
names(Earthquake)
earthquake <- cbind(Earthquake$longitude, Earthquake$latitude) 
PC <- PrincipalCircle(earthquake)
result <- GenerateCircle(PC[1:2], PC[3])
## plot
sphereplot::rgl.sphgrid(col.long = "black", col.lat = "black")
sphereplot::rgl.sphpoints(earthquake, radius = 1, col = "blue", size = 12)
sphereplot::rgl.sphpoints(result, radius = 1, col = "red", size = 6)
}
}
\keyword{ ~ exact principal circle }
