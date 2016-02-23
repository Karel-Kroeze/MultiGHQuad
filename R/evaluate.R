#' Evaluation of multivariate normal distributed expectations
#' 
#' Evaluates the posterior expectation of a (set of) parameters with a given likelihood function and multivariate normal prior distribution by Gauss-Hermite quadrature.
#'
#' The evaluated function is assumed to have a multivariate normal prior distribution, with a mean vector and covariance matrix specified in \code{\link{init.quad}}. 
#' The log-likelihood function defaults to an identity function \code{FUN(x) = 1}, which reduces the distribution under evaluation to specified prior distribution.
#' 
#' The integral under evaluation is;
#' \deqn{\int_{-\infty }^{\infty} g(X | \mu, \Sigma) \times f(X) \times X dX}{Int g(X | mu, Sigma) * f(X) * X dX}
#' where \eqn{g(X | \mu, \Sigma)}{g(X | mu, Sigma)} is the prior likelihood of X, and \eqn{f(X)} is the likelihood of X in function \code{FUN}. 
#' If left default, the result is the expectation \eqn{E(X)}, where \eqn{X \~ N(\mu, \Sigma)}{X ~ N(mu, Sigma)}.
#' 
#' 
#' Note: FUN is evaluated in a loop, vectorisation is a future possibility. FUN must return a single scalar on the natural log-scale.
#' @param FUN log likelihood function of the parameters to be estimated. 
#'     Defaults to \code{function(x) 1}, in which case only the prior likelihood is evaluated.
#' @param X Matrix of quadrature points, see \code{\link{init.quad}}. Alternatively, the list of quadrature points and weights produced by \code{\link{init.quad}}.
#' @param ... Additional arguments passed on to FUN.
#' @param W Vector of weights, or \code{NULL} (the default) if provided by \code{X}.
#' @param forcePD Logical, should the returned estimate be forced to the nearest positive definite matrix - if not already PD? If TRUE (Default: TRUE), \code{\link[Matrix]{nearPD}} is used to arrive at the closest PD matrix.
#' @param debug Logical, should we return the results of FUN?
#' @return A vector with the evaluated integrals, with attribute \code{variance} containing the (co)variance (matrix) of the estimate(s), or the positive definite matrix closest to the estiamted covariance matrix.
#' @seealso \code{\link{init.quad}} for creating quadrature points.
#' @export
#' @importFrom Matrix nearPD
#' @examples
#' ### Basic example; E(X), X ~ N(0,1)
#' grid <- init.quad(Q = 1, prior = list(mu = 0, Sigma = diag(1)))
#' eval.quad(X = grid)
#' 
#' ### Example; Rasch model person parameter
#' # E(theta), theta ~ N(0,1) * P(X = 1 | theta, beta), P is simplified rasch model
#' # set up rasch model with fixed beta, returns LL
#' rasch <- function(theta, beta, responses){
#'   p <- exp(theta - beta)/(1 + exp(theta - beta))
#'   q <- 1 - p
#'   return(log(p) * sum(responses == 1) + log(q) * sum(responses == 0))
#' }
#' 
#' # when theta == beta, P(X = 1) = .5, generate some bernoulli trials with p = .5
#' responses <- rbinom(5, 1, .5)
#' 
#' # get EAP estimate for theta, prior N(0,1)
#' eval.quad(rasch, grid, beta = 0, responses = responses)
#' 
#' # with more data, the estimate becomes more accurate, and variance decreases
#' eval.quad(rasch, grid, beta = 0, responses = rbinom(20, 1, .5))
#' eval.quad(rasch, grid, beta = 0, responses = rbinom(50, 1, .5))
#' eval.quad(rasch, grid, beta = 0, responses = rbinom(100, 1, .5))
#' 
#' ### problem; the result starts to 'snap' to the closest quadrature point when
#' # the posterior distribution is too dissimilar to the prior.
#' evals <- eval.quad(rasch, grid, beta = 0, responses = rbinom(100, 1, .5), debug = TRUE)
#' evals.values <- attr(evals, "values")
#' 
#' # posterior density after 40 items
#' p <- plot(function(x) exp(dnorm(x, log = TRUE) + 
#'                             rasch(x, beta = 1, responses = rbinom(100, 1, .5))),
#'           from = -3, to = 3)
#' 
#' # quadrature points used
#' points(grid$X, exp(grid$W)*max(p$y), pch = 20)
#' 
#' # the evaluation relies almost completely on one quadrature point, 
#' # which causes results to 'snap' to that point.
#' # we could add more quadrature points...
#' grid2 <- init.quad(Q = 1, ip = 20)
#' points(grid2$X, exp(grid2$W)*max(p$y), pch = 20, col = "grey")
#' 
#' # but if the posterior is not centered on the prior, this quickly fails:
#' p <- plot(function(x) exp(dnorm(x, log = TRUE) + 
#'                             rasch(x, beta = 2, responses = rbinom(100, 1, .5))),
#'           from = -3, to = 3)
#' points(grid2$X, exp(grid2$W)*max(p$y), pch = 20, col = "grey")
#' 
#' # additionally, adding extra quadrature points in a multidimensional 
#' # problem quickly grows out of control.
#' 
#' ### a better solution; adaptive quadrature grid.
#' # say we have an idea of where our parameter is located, through another estimator, 
#' # or a previous estimate.
#' # we can then use this to adapt where our quadrature grid should be.
#' # get an estimate;
#' responses <- rbinom(10, 1, .5)
#' est <- eval.quad(rasch, grid, beta = 2, responses = responses)
#' print( est )
#' 
#' # adapt the grid;
#' grid3 <- init.quad(Q = 1, adapt = est)
#' 
#' # grid is now much closer to posterior 
#' p <- plot(function(x) exp(dnorm(x, log = TRUE) + 
#'                             rasch(x, beta = 2, responses = rep(c(0,1), each = 20))),
#'           from = -3, to = 3)
#' points(grid3$X, exp(grid3$W)*max(p$y), pch = 20, col = "grey")
#' est <- eval.quad(rasch, grid3, beta = 2, responses = responses)
#' print(est)
eval.quad <- function(FUN = function(x) 1,
                      X = NULL,
                      ...,
                      W = NULL,
                      forcePD = TRUE,
                      debug = FALSE){
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
  
  # substitute nearest positive definite matrix
  # thanks Alwin Stegeman
  variance <- as.matrix(nearPD( variance )$mat)
  
    attr(estimate, "variance") <- variance
  
  # debug stuff
  if (debug) {
    attr(estimate, "values") <- f / p1
  }
  
  return(estimate)
}
