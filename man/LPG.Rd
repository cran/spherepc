\name{LPG}
\alias{LPG}
\title{Local principal geodesics}
\description{
Locally definded principal geodesic analysis.
}
\usage{
LPG(data, scale = 0.04, tau = scale/3, nu = 0, maxpt = 500,
seed = NULL, kernel = "indicator", thres = 1e-4, col1 = "blue", 
col2 = "green", col3 = "red")
}
\arguments{
  \item{data}{matrix or data frame consisting of spatial locations with two columns.
  Each row represents longitude and latitude.}
  \item{scale}{scale parameter for this function. The argument is the degree to which \code{LPG} expresses data locally; thus, as \code{scale} grows as the result of \code{LPG} become similar to that of \code{\link{PGA}}.}
  \item{tau}{forwarding or backwarding distance of each step. It is empirically recommended to choose a third of \code{scale}, which is the default of this argument.}
  \item{nu}{parameter to alleviate the bias of resulting curves. \code{nu} represents the viscosity of the given data and it should be selected in [0, 1). When the \code{nu} is close to 1, the curves usually swirl around, similar to the motion of a large viscous fluid. The swirling can be controlled by the argument \code{maxpt}.}
  \item{maxpt}{maximum number of points that each curve has.}
  \item{seed}{random seed number.}
  \item{kernel}{kind of kernel function. The default is the indicator kernel and alternatives are quartic or Gaussian.}
  \item{thres}{threshold of the stopping condition for the \code{IntrinsicMean} contained in the \code{LPG} function.}
  \item{col1}{color of data.}
  \item{col2}{color of points in the resulting curves.}
  \item{col3}{color of the connecting line between points in the resulting curves.}
}
\details{
Locally definded principal geodesic analysis. The result is sensitive to \code{scale} and \code{nu}, especially \code{scale} should be carefully chosen according to the structure of the given data.
}
\value{
plot and a list consisting of
\item{prin.curves}{spatial locations of points in the resulting curves.}
\item{line}{connecting line between points of \code{prin.curves}.}
\item{num.curve}{the number of the resulting curves.}
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
\donttest{
#### example 1: spiral data
n <- 500                         # the number of samples.
sigma <- 0.2                     # noise level.             
r <- runif(n) * 42                             
theta <- -pi/10 * r + 2 + sigma * rnorm(n)
x <- r * cos(theta)
y <- r * sin(theta)
simul.spiral <- cbind(x - 90, y)
LPG(simul.spiral, scale = 0.07, nu = 0.1)

#### example 2: zigzag data set
n <- 50                          # the number of samples is 6 * n = 300.
sigma <- 2                       # noise level.                 
x1 <- runif(n) * 20 - 20
y1 <- x1 + 20 + sigma * rnorm(n)
x2 <- runif(n) * 20 - 20
y2 <- - x2 + 20 + sigma * rnorm(n)
x3 <- runif(n) * 20 - 20
y3 <- x3 + 60 + sigma * rnorm(n)
x4 <- runif(n) * 20 - 20
y4 <- - x4 - 20 + sigma * rnorm(n)
x5 <- runif(n) * 20 - 20
y5 <- x5 - 20 + sigma * rnorm(n)
x6 <- runif(n) * 20 - 20
y6 <- - x6 - 60 + sigma * rnorm(n)
x <- c(x1, x2, x3, x4, x5, x6)
y <- c(y1, y2, y3, y4, y5, y6)
simul.zigzag <- cbind(x, y)
LPG(simul.zigzag, scale = 0.1, nu = 0.1)

#### example 3: earthquake data
data(Earthquake)
names(Earthquake)
earthquake <- cbind(Earthquake$longitude, Earthquake$latitude)  
LPG(earthquake, scale = 0.5, nu = 0.2, maxpt = 20)
LPG(earthquake, scale = 0.4, nu = 0.3)

#### example 4: tree data
## stem
set.seed(7)
n <- 300                  # the number of samples in stem.
sigma <- 0.1              # noise level of stem.  
lat <- 70 * runif(n) - 10
lon <- 0 + sigma * rnorm(n)
stem <- cbind(lon, lat)
## branch
n2 <- 200                 # the number of samples of each branch.
sigma2 <- 0.05            # noise level of branch. 
lon <- -20 * runif(n2)
lat <-  -lon + 10 + sigma2 * rnorm(n2)  
branch1 <- cbind(lon, lat)
lon <- 20 * runif(n2)
lat <- lon + sigma2 * rnorm(n2)  
branch3 <- cbind(lon, lat)
lon <- 20 * runif(n2)
lat <- lon + 20 + sigma2 * rnorm(n2)
branch2 <- cbind(lon, lat)
lon <- 20 * runif(n2)
lat <- lon + 40 + sigma2 * rnorm(n2)
branch5 <- cbind(lon, lat)
lon <- -20 * runif(n2)
lat <- -lon + 30 + sigma2 * rnorm(n2)
branch4 <- cbind(lon, lat)
branch <- rbind(branch1, branch2, branch3, branch4, branch5)
## sub-branches
n3 <- 20                  # the number of samples of each sub-branch.
sigma3 <- 0.01            # noise level of sub-branch.
l <- 1    # length of sub-branches
lon <- l * runif(n3) - 10
lat <- lon + 30 + sigma3 * rnorm(n3)
branches1 <- cbind(lon, lat)
lon <- -l * runif(n3) + 10
lat <- -lon + 20 + sigma3 * rnorm(n3)
branches2 <- cbind(lon, lat)
lon <- -l * runif(n3) + 10
lat <- -lon + 40 + sigma3 * rnorm(n3)
branches17 <- cbind(lon, lat)
lon <- l * runif(n3) - 14
lat <- lon + 38 + sigma3 * rnorm(n3)
branches3 <- cbind(lon, lat)
lon <- -l * runif(n3) + 14
lat <- -lon + 28 + sigma3 * rnorm(n3)
branches4 <- cbind(lon, lat)
lon <- -l * runif(n3) + 14
lat <- -lon + 48 + sigma3 * rnorm(n3)
branches18 <- cbind(lon, lat)
lon <- -l * runif(n3) - 12  
lat <- lon + 34 + sigma3 * rnorm(n3)
branches5 <- cbind(lon, lat)
lon <- l * runif(n3) + 12
lat <- -lon + 24 + sigma3 * rnorm(n3)
branches6 <- cbind(lon, lat)
lon <- l * runif(n3) + 12
lat <- -lon + 44 + sigma3 * rnorm(n3)
branches20 <- cbind(lon, lat)
lon <- -l * runif(n3) - 16
lat <- lon + 42 + sigma3 * rnorm(n3)
branches7 <- cbind(lon, lat)
lon <- l * runif(n3) + 16
lat <- -lon + 32 + sigma3 * rnorm(n3)
branches8 <- cbind(lon, lat)
lon <- l * runif(n3) + 16
lat <- -lon + 52 + sigma3 * rnorm(n3)
branches19 <- cbind(lon, lat)
lon <- l * runif(n3) + 10
lat <- -lon + 60 + sigma3 * rnorm(n3)
branches9 <- cbind(lon, lat)
lon <- -l * runif(n3) - 10
lat <- lon + 50 + sigma3 * rnorm(n3)
branches10 <- cbind(lon, lat)
lon <- -l * runif(n3) + 12
lat <- -lon + 64 + sigma3 * rnorm(n3)
branches11 <- cbind(lon, lat)
lon <- l * runif(n3) - 12
lat <- lon + 54 + sigma3 * rnorm(n3)
branches12 <- cbind(lon, lat)
lon <-  l * runif(n3) + 14
lat <- -lon + 68 + sigma3 * rnorm(n3)
branches13 <- cbind(lon, lat)
lon <- -l * runif(n3) - 14
lat <- lon + 58 + sigma3 * rnorm(n3)
branches14 <- cbind(lon, lat)
lon <- -l * runif(n3) + 16
lat <- -lon + 72 + sigma3 * rnorm(n3)
branches15 <- cbind(lon, lat)
lon <- l * runif(n3) - 16
lat <- lon + 62 + sigma3 * rnorm(n3)
branches16 <- cbind(lon, lat)
sub.branches <- rbind(branches1, branches2, branches3,
branches4,branches5, branches6, branches7, branches8, 
branches9, branches10, branches11, branches12, 
branches13, branches14, branches15, branches16, 
branches17, branches18, branches19, branches20)
tree <- rbind(stem, branch, sub.branches)
LPG(tree, scale = 0.03, nu = 0.1, seed = 7)
}
}
\keyword{ ~principal geodesic analysis }
