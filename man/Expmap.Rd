\name{Expmap}
\alias{Expmap}
\title{Exponential map}
\description{
This function performs the exponential map at (0, 0, 1) on the unit 2-sphere.
}
\usage{
Expmap(vec)
}
\arguments{
  \item{vec}{element of two-dimensional Euclidean vector space.}
}
\details{
This function performs exponential map at (0, 0, 1) on the unit sphere.
 \code{vec} is an element of the tangent plane of the unit sphere at (0, 0, 1), and the result is an element of the unit sphere in three-dimensional Euclidean space.
}
\value{
three-dimensional vector.
}
\references{
Fletcher, P. T., Lu, C., Pizer, S. M. and Joshi, S. (2004). Principal geodesic analysis for the study of nonlinear statistics of shape. IEEE Transactions on Medical Imaging, 23, 995-1005.
}
\author{
Jongmin Lee
}
\seealso{
\code{\link{Logmap}}.
}
\examples{
Expmap(c(1, 2))
}