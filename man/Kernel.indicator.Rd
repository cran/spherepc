\name{Kernel.indicator}
\alias{Kernel.indicator}
\title{Indicator kernel function}
\description{
This function returns the value of an indicator kernel function.
}
\usage{
Kernel.indicator(vec)
}
\arguments{
 \item{vec}{any length of vector.}
}
\details{
This function returns the value of an indicator kernel function. The value of kernel represents similarity from origin. The function returns a vector whose length is same as \code{vec}.
}
\value{
vector.
}
\author{
Jongmin Lee
}
\seealso{
\code{\link{Kernel.Gaussian}}, \code{\link{Kernel.quartic}}.
}
\examples{
Kernel.indicator(c(0, 1/2, 1, 2))
}

