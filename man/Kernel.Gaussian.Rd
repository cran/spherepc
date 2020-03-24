\name{Kernel.Gaussian}
\alias{Kernel.Gaussian}
\title{Gaussian kernel function}
\description{
This function returns the value of a Gaussian kernel function.
}
\usage{
Kernel.Gaussian(vec)
}
\arguments{
  \item{vec}{any length of vector.}
}
\details{
This function returns the value of a Gaussian kernel function. The value of kernel represents the similarity from origin. The function returns a vector whose length is same as \code{vec}.
}
\value{
vector.
}
\author{
Jongmin Lee 
}
\seealso{
\code{\link{Kernel.indicator}}, \code{\link{Kernel.quartic}}.
}
\examples{
Kernel.Gaussian(c(0, 1/2, 1, 2))
}

