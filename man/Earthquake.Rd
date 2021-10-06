\name{Earthquake}
\docType{data}
\alias{Earthquake}
\title{Earthquake}
\description{
It is an earthquake data from the U.S Geological Survey that collect significant earthquakes (8+ Mb magnitude) around the Pacific Ocean since the year 1900. The data are available from (\url{https://www.usgs.gov}).
Additionally, note that distribution of the data has the following features: 1) scattered, 2) curvilinear one-dimensional structure on the sphere.
}
\usage{
data(Earthquake)
}
\format{
A data frame consisting of time, latitude, longitude, depth, magnitude, etc.
}
\examples{
data(Earthquake)
names(Earthquake)
# collect spatial locations (longitude/latitude) of data.
earthquake <- cbind(Earthquake$longitude, Earthquake$latitude)  
library(rgl)
library(sphereplot)
library(geosphere)
#### example 1: principal geodesic analysis (PGA)
PGA(earthquake)
\donttest{
#### example 2: principal circle
circle <- PrincipalCircle(earthquake)          # get center and radius of principal circle
PC <- GenerateCircle(circle[1:2], circle[3])   # generate Principal circle
# plot
sphereplot::rgl.sphgrid()                      
sphereplot::rgl.sphpoints(earthquake, radius = 1, col = "blue", size = 12)
sphereplot::rgl.sphpoints(PC, radius = 1, col = "red", size = 9)

#### example 3: spherical principal curves (SPC, SPC.Hauberg)
SPC(earthquake)         # spherical principal curves.
SPC.Hauberg(earthquake) # principal curves by Hauberg on sphere.

#### example 4: local principal geodesics (LPG)
LPG(earthquake, scale = 0.5, nu = 0.2, maxpt = 20)
LPG(earthquake, scale = 0.4, nu = 0.3)  
}
}
\keyword{datasets}


 
