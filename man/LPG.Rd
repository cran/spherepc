\name{LPG}
\alias{LPG}
\title{Local principal geodesics}
\description{
Locally definded principal geodesic analysis.
}
\usage{
LPG(data, scale = 0.04, tau = scale/3, nu = 0, maxpt = 500,
seed = NULL, kernel = "indicator", thres = 1e-4, 
col1 = "blue", col2 = "green", col3 = "red")
}
\arguments{
  \item{data}{matrix or data frame consisting of spatial locations with two columns.
  Each row represents longitude and latitude (denoted by degrees).}
  \item{scale}{scale parameter for this function. The argument is the degree to which \code{LPG} expresses data locally; thus, as \code{scale} grows as the result of \code{LPG} become similar to that of the \code{\link{PGA}} function. The default is 0.4.}
  \item{tau}{forwarding or backwarding distance of each step. It is empirically recommended to choose a third of \code{scale}, which is the default of this argument.}
  \item{nu}{parameter to alleviate the bias of resulting curves. \code{nu} represents the viscosity of the given data and it should be selected in [0, 1). The default is zero. When the \code{nu} is close to 1, the curves usually swirl around, similar to the motion of a large viscous fluid. The swirling can be controlled by the argument \code{maxpt}.}
  \item{maxpt}{maximum number of points that each curve has. The default is 500.}
  \item{seed}{random seed number.}
  \item{kernel}{kind of kernel function. The default is the indicator kernel and alternatives are quartic or Gaussian.}
  \item{thres}{threshold of the stopping condition for the \code{IntrinsicMean} contained in the \code{LPG} function. The default is 1e-4.}
  \item{col1}{color of data. The default is blue.}
  \item{col2}{color of points in the resulting principal curves. The default is green.}
  \item{col3}{color of the resulting curves. The default is red.}
}
\details{
Locally definded principal geodesic analysis. The result is sensitive to \code{scale} and \code{nu}, especially \code{scale} should be carefully chosen according to the structure of the given data.
}
\value{
plot and a list consisting of
\item{prin.curves}{spatial locations (represented by degrees) of points in the resulting curves.}
\item{line}{connecting lines between points in \code{prin.curves}.}
\item{num.curves}{the number of the resulting curves.}
}
\author{
Jongmin Lee
}
\seealso{
\code{\link{PGA}}, \code{\link{SPC}}, \code{\link{SPC.Hauberg}}.
}
\examples{
library(rgl)
library(sphereplot)
library(geosphere)

#### example 1: spiral data
## longitude and latitude are expressed in degrees
set.seed(40)
n <- 900                                    # the number of samples
sigma1 <- 1; sigma2 <- 2.5;                 # noise levels
radius <- 73; slope <- pi/16                # radius and slope of spiral
## polar coordinate of (longitude, latitude)
r <- runif(n)^(2/3) * radius; theta <- -slope * r + 3 
## transform to (longitude, latitude)
correction <- (0.5 * r/radius + 0.3)        # correction of noise level
lon1 <- r * cos(theta) + correction * sigma1 * rnorm(n)
lat1 <- r * sin(theta) + correction * sigma1 * rnorm(n)
lon2 <- r * cos(theta) + correction * sigma2 * rnorm(n)
lat2 <- r * sin(theta) + correction * sigma2 * rnorm(n)
spiral1 <- cbind(lon1, lat1); spiral2 <- cbind(lon2, lat2)
## plot spiral data
rgl.sphgrid(col.lat = 'black', col.long = 'black')
rgl.sphpoints(spiral1, radius = 1, col = 'blue', size = 12)
## implement the LPG to (noisy) spiral data
\donttest{LPG(spiral1, scale = 0.06, nu = 0.1, seed = 100)}
\donttest{LPG(spiral2, scale = 0.12, nu = 0.1, seed = 100)}
#### example 2: zigzag data set
set.seed(10)
n <- 50                                # the number of samples is 6 * n = 300
sigma1 <- 2; sigma2 <- 5               # noise levels                   
x1 <- x2 <- x3 <- x4 <- x5 <- x6 <- runif(n) * 20 - 20
y1 <- x1 + 20 + sigma1 * rnorm(n); y2 <- -x2 + 20 + sigma1 * rnorm(n)
y3 <- x3 + 60 + sigma1 * rnorm(n); y4 <- -x4 - 20 + sigma1 * rnorm(n)
y5 <- x5 - 20 + sigma1 * rnorm(n); y6 <- -x6 - 60 + sigma1 * rnorm(n)
x <- c(x1, x2, x3, x4, x5, x6); y <- c(y1, y2, y3, y4, y5, y6)
simul.zigzag1 <- cbind(x, y)
## plot zigzag data
sphereplot::rgl.sphgrid(col.lat = 'black', col.long = 'black')
sphereplot::rgl.sphpoints(simul.zigzag1, radius = 1, col = 'blue', size = 12)
## implement the LPG to zigzag data
\donttest{LPG(simul.zigzag1, scale = 0.1, nu = 0.1, maxpt = 45, seed = 50)}

## noisy zigzag data
set.seed(10)
z1 <- z2 <- z3 <- z4 <- z5 <- z6 <- runif(n) * 20 - 20
w1 <- z1 + 20 + sigma2 * rnorm(n); w2 <- -z2 + 20 + sigma2 * rnorm(n)
w3 <- z3 + 60 + sigma2 * rnorm(n); w4 <- -z4 - 20 + sigma2 * rnorm(n)
w5 <- z5 - 20 + sigma2 * rnorm(n); w6 <- -z6 - 60 + sigma2 * rnorm(n)
z <- c(z1, z2, z3, z4, z5, z6); w <- c(w1, w2, w3, w4, w5, w6)
simul.zigzag2 <- cbind(z, w)
## implement the LPG to noisy zigzag data
\donttest{LPG(simul.zigzag2, scale = 0.2, nu = 0.1, maxpt = 18, seed = 20)}


#### example 3: Doubly circular data set
set.seed(30)
n <- 200
sigma <- 1
x1 <- 40 * runif(n) - 20
y1 <- (-x1^2 + 400)^(1/2) + 30 + sigma * rnorm(n)
x2 <- 40 * runif(n) - 20
y2 <- -(-x2^2 + 400)^(1/2) + 30 + sigma * rnorm(n)
y3 <- 40 * runif(n) + 10
x3 <- -(-y3^2 + 60 * y3 - 500)^(1/2) + sigma * rnorm(n)
y4 <- 40 * runif(n) + 10
x4 <- (-y4^2 + 60 * y4 - 500)^(1/2) + sigma * rnorm(n)
Dc1 <- cbind(c(x1, x2, x3, x4), c(y1, y2, y3, y4))
z1 <- 40 * runif(n) - 20
w1 <- (400 - z1^2)^(1/2) - 20 + sigma * rnorm(n)
z2 <- 40 * runif(n) - 20
w2 <- -(400 - z2^2)^(1/2) - 20 + sigma * rnorm(n)
w3 <- -40 * runif(n)
z3 <- (-w3^2 - 40 * w3)^(1/2) + sigma * rnorm(n)
w4 <- -40 * runif(n)
z4 <- -(-w4^2 - 40 * w4)^(1/2) + sigma * rnorm(n)
Dc2 <- cbind(c(z1, z2, z3, z4), c(w1, w2, w3, w4))
Dc <- rbind(Dc1, Dc2)
\donttest{LPG(Dc, scale = 0.15, nu = 0.1, maxpt = 22,)}


#### example 4: real earthquake data
data(Earthquake)
names(Earthquake)
earthquake <- cbind(Earthquake$longitude, Earthquake$latitude)  
LPG(earthquake, scale = 0.5, nu = 0.2, maxpt = 20)
LPG(earthquake, scale = 0.4, nu = 0.3)


#### example 5: tree data
## tree consists of stem, branches and subbranches

## generate stem
set.seed(10)
n1 <- 200; n2 <- 100; n3 <- 15 # the number of samples in stem, a branch, and a subbranch
sigma1 <- 0.1; sigma2 <- 0.05; sigma3 <- 0.01 # noise levels
noise1 <- sigma1 * rnorm(n1); noise2 <- sigma2 * rnorm(n2); noise3 <- sigma3 * rnorm(n3)
l1 <- 70; l2 <- 20; l3 <- 1            # length of stem, branches, and subbranches
rep1 <- l1 * runif(n1)                 # repeated part of stem
stem <- cbind(0 + noise1, rep1 - 10)

## generate branch
rep2 <- l2 * runif(n2)                 # repeated part of branch
branch1 <- cbind(-rep2, rep2 + 10 + noise2); branch2 <- cbind(rep2, rep2 + noise2)
branch3 <- cbind(rep2, rep2 + 20 + noise2); branch4 <- cbind(rep2, rep2 + 40 + noise2)
branch5 <- cbind(-rep2, rep2 + 30 + noise2)
branch <- rbind(branch1, branch2, branch3, branch4, branch5)

## generate subbranches
rep3 <- l3 * runif(n3)                 # repeated part in subbranches
branches1 <- cbind(rep3 - 10, rep3 + 20 + noise3)
branches2 <- cbind(-rep3 + 10, rep3 + 10 + noise3)
branches3 <- cbind(rep3 - 14, rep3 + 24 + noise3)
branches4 <- cbind(-rep3 + 14, rep3 + 14 + noise3)
branches5 <- cbind(-rep3 - 12, -rep3 + 22 + noise3)
branches6 <- cbind(rep3 + 12, -rep3 + 12 + noise3)
branches7 <- cbind(-rep3 - 16, -rep3 + 26 + noise3)
branches8 <- cbind(rep3 + 16, -rep3 + 16 + noise3)
branches9 <- cbind(rep3 + 10, -rep3 + 50 + noise3)
branches10 <- cbind(-rep3 - 10, -rep3 + 40 + noise3)
branches11 <- cbind(-rep3 + 12, rep3 + 52 + noise3)
branches12 <- cbind(rep3 - 12, rep3 + 42 + noise3)
branches13 <- cbind(rep3 + 14, -rep3 + 54 + noise3)
branches14 <- cbind(-rep3 - 14, -rep3 + 44 + noise3)
branches15 <- cbind(-rep3 + 16, rep3 + 56 + noise3)
branches16 <- cbind(rep3 - 16, rep3 + 46 + noise3)
branches17 <- cbind(-rep3 + 10, rep3 + 30 + noise3)
branches18 <- cbind(-rep3 + 14, rep3 + 34 + noise3)
branches19 <- cbind(rep3 + 16, -rep3 + 36 + noise3)
branches20 <- cbind(rep3 + 12, -rep3 + 32 + noise3)
sub.branches <- rbind(branches1, branches2, branches3, branches4, branches5, branches6,
+    branches7, branches8, branches9, branches10, branches11, branches12, branches13,
+    branches14, branches15, branches16, branches17, branches18, branches19, branches20)

## tree consists of stem, branch and subbranches
tree <- rbind(stem, branch, sub.branches)

## plot tree data
sphereplot::rgl.sphgrid(col.lat = 'black', col.long = 'black')
sphereplot::rgl.sphpoints(tree, radius = 1, col = 'blue', size = 12)

## implement the LPG function to tree data
\donttest{LPG(tree, scale = 0.03, nu = 0.2, seed = 10)}
}
\keyword{ ~principal geodesic analysis }
\keyword{ ~principal curves }