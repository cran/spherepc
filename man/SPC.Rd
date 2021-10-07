\name{SPC}
\alias{SPC}
\title{Spherical principal curves}
\description{
This function fits a spherical principal curve.
}
\usage{
SPC(data, q = 0.1, T = nrow(data), step.size = 1e-3, maxit = 10, 
type = "Intrinsic", thres = 0.1, deletePoints = FALSE, plot.proj = FALSE, 
kernel = "quartic", col1 = "blue", col2 = "green", col3 = "red")
}
\arguments{
  \item{data}{matrix or data frame consisting of spatial locations with two columns.
  Each row represents a longitude and latitude (denoted by degrees).}
  \item{q}{numeric value of the smoothing parameter. The parameter plays the same role, as the bandwidth does in kernel regression, in the \code{SPC} function. The value should be a numeric value between 0.01 and 0.5. The default is 0.1.}
  \item{T}{the number of points making up the resulting curve. The default is 1000.}
  \item{step.size}{step size of the \code{PrincipalCircle} function. The default is 0.001. The resulting principal circle is used by an initial estimate of the \code{SPC}.}
  \item{maxit}{maximum number of iterations. The default is 30.}
  \item{type}{type of mean on the sphere. The default is "Intrinsic" and the other choice is "Extrinsic".}
  \item{thres}{threshold of the stopping condition. The default is 0.1}
  \item{deletePoints}{logical value. The argument is an option of whether to delete points or not. If \code{deletePoints} is FALSE, this function leaves the points in curves which do not have adjacent data for each expectation step. As a result, the function usually returns a closed curve, i.e. a curve without endpoints. If \code{deletePoints} is TRUE, this function deletes the points in curves which do not have adjacent data for each expectation step. As a result, The \code{SPC} function usually returns an open curve, i.e. a curve with endpoints. The default is FALSE.}
  \item{plot.proj}{logical value. If the argument is TRUE, the projection line for each data is plotted. The default is FALSE.}
  \item{kernel}{kind of kernel function. The default is quartic kernel and alternatives are indicator or Gaussian.}
  \item{col1}{color of data. The default is blue.}
  \item{col2}{color of points in the principal curves. The default is green.}
  \item{col3}{color of resulting principal curves. The default is red.}
}
\details{
This function fits a spherical principal curves, and requires to load the 'rgl', 'sphereplot', and 'geosphere' R packages.
}
\value{
plot and a list consisting of
  \item{prin.curves}{spatial locations (denoted by degrees) of points in the resulting principal curves.}
  \item{line}{connecting line bewteen points of \code{prin.curves}.}
  \item{converged}{whether or not the algorithm converged.}
  \item{iteration}{the number of iterations of the algorithm.}
  \item{recon.error}{sum of squared distances from the data to their projections.}
  \item{num.dist.pt}{the number of distinct projections.}
}
\references{
Jang-Hyun Kim, Jongmin Lee and Hee-Seok Oh. (2020). Spherical Principal Curves <arXiv:2003.02578>.

Jongmin Lee, Jang-Hyun Kim and Hee-Seok Oh. (2021). Spherical principal curves. IEEE Transactions on Pattern Analysis and Machine Intelligence, 43, 2165-2171. <doi.org/10.1109/TPAMI.2020.3025327>.
}
\author{
Jongmin Lee
}
\note{
This function requires to load 'rgl', 'sphereplot', and 'geosphere' R packages.
}
\seealso{
\code{\link{SPC.Hauberg}}.
}
\examples{
library(rgl)
library(sphereplot)
library(geosphere)

#### example 2: waveform data
n <- 200
alpha <- 1/3; freq <- 4                           # amplitude and frequency of wave
sigma1 <- 2; sigma2 <- 10                         # noise levels  
lon <- seq(-180, 180, length.out = n)             # uniformly sampled longitude
lat <- alpha * 180/pi * sin(freq * lon * pi/180) + 10.        # latitude vector
## add Gaussian noises on latitude vector
lat1 <- lat + sigma1 * rnorm(length(lon)); lat2 <- lat + sigma2 * rnorm(length(lon))
wave1 <- cbind(lon, lat1); wave2 <- cbind(lon, lat2)
## implement SPC to the (noisy) waveform data
\donttest{SPC(wave1, q = 0.05)}
\donttest{SPC(wave2, q = 0.05)}

#### example 1: earthquake data
data(Earthquake)
names(Earthquake)
earthquake <- cbind(Earthquake$longitude, Earthquake$latitude)   
\donttest{SPC(earthquake, q = 0.1)
## options 1: plot the projection lines (use option of plot.proj = TRUE)
SPC(earthquake, q = 0.1, plot.proj = TRUE)
## option 2: open principal curves (use option of deletePoints = TRUE)
SPC(earthquake, q = 0.04, deletePoints = TRUE)}

}
\keyword{ ~principal curves }

