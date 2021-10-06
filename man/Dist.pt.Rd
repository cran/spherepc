\name{Dist.pt}
\alias{Dist.pt}
\title{The number of distinct points.}
\description{
This function calculates the number of distinct point in the given data.
}
\usage{
Dist.pt(data)
}
\arguments{
  \item{data}{matrix or dataframe consisting of spatial locations with two columns.
  Each row represents longitude and latitude.}
}
\details{
This function calculates the number of distinct points in the given data.
}
\value{
a numeric.
}
\author{
Jongmin Lee
}
\examples{
Dist.pt(rbind(c(0, 0), c(0, 1), c(1, 0), c(1, 1), c(0, 0)))
}
