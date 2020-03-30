\name{SPC.Hauberg}
\alias{SPC.Hauberg}
\title{principal curves by Hauberg on the sphere}
\description{
This function fits a principal curve by Hauberg on the sphere.
}
\usage{
SPC.Hauberg(data, q = 0.1, T = nrow(data), step.size = 1e-3, maxit = 30, 
type = "Intrinsic", thres = 1e-2, deletePoints = FALSE, plot.proj = FALSE, 
kernel = "quartic", col = c("blue", "green", "red", "black"), size = c(12, 6, 6))
}
\arguments{
  \item{data}{matrix or data frame consisting of spatial locations with two columns. Each row represents longitude and latitude.}
  \item{q}{numeric value of the smoothing parameter. Intuitively speaking, the role of this argument is similar to the that of bandwidth for kernel regression. The value should be a numeric value between 0.01 and 0.5. The default is 0.1.}
  \item{T}{the number of points in the resulting curve.}
  \item{step.size}{step size of the \code{PrincipalCircle}. The resulting principal circle is used by an initialization of the \code{SPC}.}
  \item{maxit}{maximum number of iterations.}
  \item{type}{type of mean on sphere. The default is "Intrinsic" and the alternative is "extrinsic".}
  \item{thres}{threshold of the stopping condition.}
  \item{deletePoints}{logical value. The argument is an option of whether to delete points or not. If \code{deletePoints} is FALSE, this function leaves the points in curves which do not have adjacent data for each expectation step. As a result, the function usually returns a closed curve, i.e., a curve without endpoints. If \code{deletePoints} is TRUE, this function deletes the points in curves which do not have adjacent data for each expectation step. As a result, The \code{SPC} function usually returns an open curve, i.e., a curve with endpoints. The default is FALSE.}
  \item{plot.proj}{logical value. If the argument is TRUE, the projection line for each data is plotted. The default is FALSE.}
  \item{kernel}{kind of kernel function. The default is quartic kernel and alternatives are indicator or Gaussian.}
  \item{col}{four-dimensional vector which represents colors of data, points in the resulting curves, the connecting line between points in the resulting curves, and projection lines, respectively.}
  \item{size}{three-dimensional vector which represents sizes of data, points in the principal curves, and the connecting line between points in the curves, respectively.}
}
\details{
This function fits a principal curve proposed by Hauberg on the sphere, and requires to load the 'rgl', 'sphereplot', and 'geosphere' R packages.
}
\value{
plot and a list consisting of
  \item{prin.curves}{spatial points of in the resulting principal curves.}
  \item{line}{connecting line bewteen points of \code{prin.curves}.}
  \item{converged}{whether or not the algorithm converged.}
  \item{iteration}{the number of iterations of the algorithm.}
  \item{recon.error}{sum of squared distances from the data to their projections.}
  \item{num.dist.pt}{the number of distinct projections.}
  \item{plot}{plotting of the data and principal curves.}
}
\references{
Hauberg, S. (2016). Principal curves on Riemannian manifolds. IEEE Transactions on Pattern Analysis and Machine Intelligence, 38, 1915-1921.

Jang-Hyun Kim, Jongmin Lee, and Hee-Seok Oh. (2020). Spherical Principal Curves <arXiv:2003.02578>.
}
\author{
Jongmin Lee
}
\note{
This function requires to load 'rgl', 'sphereplot', and 'geosphere' R packages.
}
\seealso{
\code{\link{SPC}}, \code{\link{Proj.Hauberg}}.
}
\examples{
library(rgl)
library(sphereplot)
library(geosphere)
#### example 1: earthquake data
data(Earthquake)
names(Earthquake)
earthquake <- cbind(Earthquake$longitude, Earthquake$latitude)
\donttest{SPC.Hauberg(earthquake, q = 0.1)}
#### example 2: waveform data
n <- 200
alpha <- 1/3    # amplitude
freq <- 4       # frequency
sigma <- 2
lon <- seq(-180, 180, length.out = n)
lat <- alpha * 180/pi * sin(lon * pi/180 * freq) + 10 + sigma * rnorm(length(lon))
wave <- cbind(lon, lat)
\donttest{SPC.Hauberg(wave, q = 0.05)}
}
\keyword{ ~principal curves }
\keyword{ ~spherical surface }
