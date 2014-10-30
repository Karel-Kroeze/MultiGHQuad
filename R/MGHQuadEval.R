#' Evaluation of multivariate normal distributed integral
#' 
#' Evaluates a given function with a (built-in) multivariate normal distribution with Gauss-Hermite quadrature.
#'
#' The evaluated function is assumed to have a multivariate normal distribution, with a given mean vector and covariance matrix. 
#' The default identity function \code{function(x) 1} reduces to an integral over a multivariate normal distribution with mean vector \code{mu} and covariance matrix \code{Sigma}.
#' 
#' @param Q Number of dimensions. Defaults to 2. Only required when \code{mu} and \code{Sigma} are not provided.
#' @param mu Mean vector, defaults to rep(0,Q), the zero vector of length Q.
#' @param Sigma Covariance matrix, defaults to diag(Q), the identity matrix of rank Q.
#' @param X Matrix of quadrature points, see \code{\link{MGHQuadPoints}}. Alternatively, the list of quadrature points and weights produced by \code{\link{MGHQuadPoints}}.
#' @param W Vector of weights, or \code{NULL} if provided by \code{X}.
#' @param ... Additional arguments passed on to FUN.
#' @return A vector with the evaluated integrals.
#' @seealso \code{\link{MGHQuadPoints}} for creating quadrature points.
#' @export
#' @examples
#' quadPoints <- MGHQuadPoints(Q=3)
#' # expected value of 3-dimensional multivariate normal distribution: N(0,1). (Since mean is currently fixed at zero, this is always zero.)
#' integral <- MGHQuadEval(Q=3,X=quadPoints) 
#' integral
#' round(integral)

MGHQuadEval <- function(FUN = function(x) 1,Q=2,mu=rep(0,Q),Sigma=diag(Q),X=NULL,W=NULL,...){
  if (is.list(X)){
    W <- X$W
    X <- X$X
  }
  if (is.null(X) | is.null(W)) stop("Quadrature points and weights are required. See MGHQuadPoints.", call.=F)

  ipq <- length(W)
  aux <- numeric(ipq)
  
  # main loop
  for (i in 1:ipq){
    aux[i] <- FUN(X[i,],...) * W[i]
  }
  
  # normalizing constant
  p1 <- sum(aux)
  
  # multiply integrals with x values, sum over columns, divide by normalizing constant.
  out <- (rep(1,ipq) %*% (aux * X)) / p1
  
  # return computed integral.
  return(out)
}


