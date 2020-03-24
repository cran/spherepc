\name{Cal.recon}
\alias{Cal.recon}
\title{Calculating reconstruction error}
\description{
This function calculates reconstruction error.
}
\usage{
Cal.recon(data, line)
}
\arguments{
  \item{data}{matrix or data frame consisting of spatial locations with two columns.
  Each row represents longitude and latitude.}
  \item{line}{longitude and latitude of a line as a matrix or data frame with two columns.}
}
\details{
This function calculates reconstruction error from the data to the line. Longitude should range from -180 to 180 and latitude from -90 to 90. This function requires to load 'geosphere' R package.
}
\value{
summation of squared distance from the data to the line on the unit sphere.
}
\author{
Jongmin Lee
}
\examples{
library(geosphere)           # This function needs to load 'geosphere' R packages.
data <- rbind(c(0, 0), c(50, -10), c(100, -70))
line <- rbind(c(30, 30), c(-20, 50), c(50, 80))
Cal.recon(data, line)
}
