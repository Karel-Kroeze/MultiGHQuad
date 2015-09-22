#' Multidimensional Gauss-Hermite Quadrature
#' 
#' Functions to perform n-dimensional numerical integration on n parameters with a multivariate normal prior distribution.
#' 
#' Use \code{\link{init.quad}} to generate a quadrature grid, and \code{\link{eval.quad}} to evaluate the integral.
#' Evaluation is performed with Gauss-Hermite quadrature, with a prior distribution that can be specified to any multivariate normal. 
#' Additionally, the grid can be adapted to any multivariate normal distribution - that is known to be close(r) 
#' to the posterior distribution under evaluation.
#'   
#' @name MultiGHQuad-package
#' @aliases MultiGHQuad, MultiGHQuad-package
#' @docType package
#' @seealso \code{\link{init.quad}}, \code{\link{eval.quad}}
#' @author Karel A Kroeze, \email{k.a.kroeze@@utwente.nl}
#' @references Jaeckel, P. (2005). \emph{A note on multivariate Gauss-Hermite quadrature}. London: ABN-Amro. Retrieved from http://www.pjaeckel.webspace.virginmedia.com/ANoteOnMultivariateGaussHermiteQuadrature.pdf
#' @references Bock, R. D., & Mislevy, R. J. (1982). Adaptive EAP Estimation of Ability in a Microcomputer Environment. \emph{Applied Psychological Measurement, 6}(4), 431-444. http://doi.org/10.1177/014662168200600405
NULL