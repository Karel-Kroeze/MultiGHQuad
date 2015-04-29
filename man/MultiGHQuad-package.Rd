\name{MultiGHQuad-package}
\alias{MultiGHQuad-package}
\alias{MultiGHQuad}
\docType{package}
\title{
Multidimensional Gauss-Hermite quadrature
}
\description{
Chained numerical integration of multivariate normal distributed parameters.
}
\details{
\tabular{ll}{
Package: \tab MultiGHQuad\cr
Type: \tab Package\cr
Version: \tab 0.2\cr
Date: \tab 2014-10-27\cr
License: \tab GPL-3\cr
}
\code{\link{init.quad}} to generate quadrature points, \code{\link{eval.quad}} to evaluate the integral. 
}
\author{
Karel Kroeze

Maintainer: Karel Kroeze <k.a.kroeze@student.utwente.nl>
}
\references{
Jaeckel, P. (2005). \emph{A note on multivariate Gauss-Hermite quadrature}. London: ABN-Amro. Retrieved from http://www.pjaeckel.webspace.virginmedia.com/ANoteOnMultivariateGaussHermiteQuadrature.pdf
}

\keyword{ package }
\seealso{
\code{\link[fastGHQuad:fastGHQuad-package]{fastGHQuad}}
}
\examples{
# Leaving FUN blank, expected value on each dimension 
# of a Q-dimensional multivariate normal distribution. 
quadPoints <- init.quad(Q=3)
eval.quad(X = quadPoints)
}
