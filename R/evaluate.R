#' Evaluation of multivariate normal distributed integral
#' 
#' Evaluates a given function with a (built-in) multivariate normal distribution with Gauss-Hermite quadrature.
#'
#' The evaluated function is assumed to have a multivariate normal distribution, with a given mean vector and covariance matrix. 
#' The default identity function \code{function(x) 1} reduces to an integral over a multivariate normal distribution with mean vector \code{mu} and covariance matrix \code{Sigma}.
#' 
#' @param FUN LOG likelihood function of the parameters to be estimated. 
#'     Defaults to \code{funtion(x) 1}, in which case only the built-in multivariate normal pdf is evaluated.
#' @param X Matrix of quadrature points, see \code{\link{init.quad}}. Alternatively, the list of quadrature points and weights produced by \code{\link{init.quad}}.
#' @param W Vector of weights, or \code{NULL} if provided by \code{X}.
#' @param debug Logical, should we return evaluated FUN results?
#' @param ... Additional arguments passed on to FUN.
#' @return A vector with the evaluated integrals, as well as a covariance matrix attribute.
#' @seealso \code{\link{init.quad}} for creating quadrature points.
#' @export
#' @examples
#' quadPoints <- init.quad(Q=3)
#' # expected value of 3-dimensional multivariate normal distribution: N(0,1). 
#' # (Since mean is currently fixed at zero, this is always zero.)
#' integral <- eval.quad(X=quadPoints) 
#' integral
#' round(integral)

eval.quad <- function(FUN = function(x) 1, X = NULL, W = NULL, debug = FALSE, ...){
  # allow list input
  if (is.list(X)){
    W <- X$W
    X <- X$X
  }
  
  # require points + weights
  if (is.null(X) | is.null(W)) stop("Quadrature points and weights are required. See init.gauss.", call.=F)
  
  # find or create function
  FUN <- match.fun(FUN)
  
  # prepare working vars
  Q <- ncol(X)
  ipq <- length(W)
  f <- numeric(ipq)
  
  # main loop
  for (i in 1:ipq){
    f[i] <- FUN(X[i,],...) + W[i]
  }
  
  # some numerical safeguards
  m <- 700 - max(f)
  f <- f + m
  
  # back to normal scale
  f <- exp(f)
  
  # normalizing constant
  p1 <- sum(f)
  
  # multiply integrals with x values, sum over columns, divide by normalizing constant.
  estimate <- colSums(f * X) / p1
  
  # variance estimate (Bock & Mislevy, 1982)
  variance <- matrix(0, Q, Q)
  for (i in 1:ipq){
    deviation <- X[i,] - estimate
    variance <- variance + ( deviation %*% t(deviation) * f[i] / p1 )
  }
  
  attr(estimate, "variance") <- variance
  
  # debug stuff
  if (debug) {
    attr(estimate, "values") <- f
  }
  
  return(estimate)
}
