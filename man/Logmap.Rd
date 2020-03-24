\name{Logmap}
\alias{Logmap}
\title{Logarithm map}
\description{
This function performs the logarithm map at (0, 0, 1) on the unit sphere.
}
\usage{
Logmap(vec)
}
\arguments{
  \item{vec}{an element of the unit sphere in three-dimensional Euclidean space.}
}
\details{
This function performs the logarithm map at (0, 0, 1) on the unit sphere. Note that, \code{vec} is normalized to be contained in the unit sphere.
}
\value{
two-dimensional vector.
}
\references{
Fletcher, P. T., Lu, C., Pizer, S. M. and Joshi, S. (2004). Principal geodesic analysis for the study of nonlinear statistics of shape. IEEE Transactions on Medical Imaging, 23, 995-1005.
}
\author{
Jongmin Lee
}
\seealso{
\code{\link{Expmap}}.
}
\examples{
Logmap(c(1/sqrt(2), 1/sqrt(2), 0))
}

