\name{Kernel.quartic}
\alias{Kernel.quartic}
\title{Quartic kernel function}
\description{
This function returns the value of a quartic kernel function.
}
\usage{
Kernel.quartic(vec)
}
\arguments{
  \item{vec}{any length of vector.}
}
\details{
This function returns the value of quartic kernel function. The value of kernel represents similarity from origin. The function returns a vector whose length is same as \code{vec}.
}
\value{
vector.
}
\author{
Jongmin Lee 
}
\seealso{
\code{\link{Kernel.Gaussian}}, \code{\link{Kernel.indicator}}.
}
\examples{
Kernel.quartic(c(0, 1/2, 1, 2))
}

