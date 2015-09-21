#' Evaluation of multivariate normal distributed expectations
#' 
#' Evaluates the expectation of a (set of) parameters with a given function under a (built-in) multivariate normal distribution by Gauss-Hermite quadrature.
#'
#' The evaluated function is assumed to have a multivariate normal distribution, with a given mean vector and covariance matrix. 
#' The default identity function \code{FUN(x) = 1} reduces a multivariate normal distribution with mean vector \code{mu} and covariance matrix \code{Sigma}.
#' 
#' The integral under evaluation is;
#' \deqn{\int_{-\Infty}^{\Infty} N(\mu, \Sigma) \times FUN(X) \times X dX}{Int Prior * FUN(X) * X dX}
#' 
#' Hence, if left default, the result is the expectation \eqn{E(N(\mu,\Sigma))}{ of N(mu, Sigma)}.
#' 
#' Note: FUN is evaluated in a loop, vectorisation is a future possibility. FUN must return a single scalar on the log-scale.
#' 
#' 
#' @param FUN log likelihood function of the parameters to be estimated. 
#'     Defaults to \code{funtion(x) 1}, in which case only the built-in multivariate normal pdf is evaluated.
#' @param X Matrix of quadrature points, see \code{\link{init.quad}}. Alternatively, the list of quadrature points and weights produced by \code{\link{init.quad}}.
#' @param W Vector of weights, or \code{NULL} if provided by \code{X}.
#' @param debug Logical, should we return evaluated FUN results?
#' @param ... Additional arguments passed on to FUN.
#' @return A vector with the evaluated integrals, as well as a covariance matrix attribute.
#' @seealso \code{\link{init.quad}} for creating quadrature points.
#' @export
eval.quad <- function(FUN = function(x) 1, X = NULL, ..., W = NULL, debug = FALSE){
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
