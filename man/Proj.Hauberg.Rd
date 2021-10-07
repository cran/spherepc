\name{Proj.Hauberg}
\alias{Proj.Hauberg}
\title{Projecting the nearest point}
\description{
This function performs the approximated projection for each data.
}
\usage{
Proj.Hauberg(data, line)
}
\arguments{
  \item{data}{matrix or data frame consisting of spatial locations with two columns.
  Each row represents longitude and latitude.}
  \item{line}{longitude and latitude of line as a matrix or data frame with two columns.}
}
\details{
This function returns the nearest points in \code{line} for each point in the data.
The function requires to load the 'geosphere' R package.
}
\value{
matrix consisting of spatial locations with two columns.
}
\references{
Hauberg, S. (2016). Principal curves on Riemannian manifolds. IEEE Transactions on Pattern Analysis and Machine Intelligence, 38, 1915-1921.

Jang-Hyun Kim, Jongmin Lee and Hee-Seok Oh. (2020). Spherical Principal Curves <arXiv:2003.02578>.

Jongmin Lee, Jang-Hyun Kim and Hee-Seok Oh. (2021). Spherical principal curves. IEEE Transactions on Pattern Analysis and Machine Intelligence, 43, 2165-2171. <doi.org/10.1109/TPAMI.2020.3025327>.
}
\author{
Jongmin Lee
}
\seealso{
\code{\link{SPC.Hauberg}}
}
\examples{
library(geosphere)
Proj.Hauberg(rbind(c(0, 0), c(10, -20)), rbind(c(50, 10), c(40, 20), c(30, 30)))
}
