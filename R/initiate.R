#' Q-dimensional grid of quadrature points.
#' 
#' Creates a flattened, rotated grid that incorporates correlation through an eigenvalue decomposition of the covariance matrix.
#'
#' Creates a Q-dimensional grid by calling \code{\link{expand.grid}} on Q vectors of unidimensional quadrature points obtained with \code{\link[fastGHQuad]{gaussHermiteData}}.
#' The grid is then corrected for a prior distribution, and can optionally be adapted around a previous estimate. The resultant grid can be pruned to remove quadrature points that are unlikely to add information.
#' 
#' @param Q Number of dimensions. Defaults to 2. Only required when \code{mu} and \code{Sigma} are not provided.
#' @param prior List of prior mean \code{mu}, = \code{vector}, and covariance matrix \code{Sigma} = \code{matrix}, defaults to zero vector and identity matrix respectively.
#' @param adapt List of adaptive mean \code{mu}, = \code{vector}, and covariance matrix \code{Sigma} = \code{matrix}, if \code{NULL} no adaptation is used. Defaults to NULL.
#' @param ip Number of quadrature points \emph{per dimension}. Defaults to 6. Note that the total number of quadrature points is \code{ip^Q}.
#' @param prune Logical, should quadrature points with a very low weight be removed? Defaults to false. See details.
#' @param forcePD Logical, should \code{adapt} and \code{prior} arguments be forced to the neares positive definite matrix - if not already PD? If TRUE (Default: FALSE), \code{\link[Matrix]{nearPD}} is used to arrive at the closest PD matrix. 
#' @param debug Logical, draws debugging plots when true.
#' @return A list with a matrix \code{X} of \code{ip^Q} by \code{Q} quadrature points and a vector \code{W} of length \code{ip^Q} associated weights.
#' @seealso \code{\link[fastGHQuad]{gaussHermiteData}}, used to create unidimensional quadrature points, and \code{\link{eval.quad}} for evaluating the integral.
#' @export
#' @importFrom mvtnorm rmvnorm
#' @importFrom Matrix nearPD
#' @examples 
#' ### basic quadrature grid /w pruning.
#' mu <- c(0,0)
#' sigma <- matrix(c(1,.5,.5,1),2,2)
#' grid <- init.quad(Q = 2, prior = list(mu = mu, Sigma = sigma), ip = 10, prune = FALSE)
#' grid2 <- init.quad(Q = 2, prior = list(mu = mu, Sigma = sigma), ip = 10, prune = TRUE)
#' library(mvtnorm)
#' normal <- rmvnorm(1000, mu, sigma)
#' # noise
#' plot(normal, xlim = c(-6,6), ylim = c(-6,6), pch = 19, col = rgb(0,0,0,.5))
#' # full quad grid
#' points(grid$X, cex = exp(grid$W)/max(exp(grid$W))*4, col = 'red', pch = 20)
#' # pruned quad grid
#' points(grid2$X, cex = exp(grid2$W)/max(exp(grid2$W))*4, col = 'green', pch = 20)
#' 
#' 
#' ### Adaptive quadrature grid
#' prior <- list(mu = c(0,0), Sigma = matrix(c(1,.5,.5,1),2,2))
#' adapt <- list(mu = c(-2,2), Sigma = prior$Sigma / 2)
#' grid <- init.quad(Q = 2, prior, ip = 10, prune = FALSE)
#' library(mvtnorm)
#' normal <- rmvnorm(1000, adapt$mu, adapt$Sigma)
#' # noise, centered at (-2, 2)
#' plot(normal, xlim = c(-6,6), ylim = c(-6,6), pch = 19, col = rgb(0,0,0,.5))
#' # initial quad grid, centered at (0, 0)
#' points(grid$X, cex = exp(grid$W)/max(exp(grid$W))*4, col = 'red', pch = 20)
#' # adapted grid
#' grid2 <- init.quad(Q =2, prior, adapt = adapt, ip = 10, prune = TRUE)
#' points(grid2$X, cex = exp(grid2$W)/max(exp(grid2$W))*4, col = 'green', pch = 20)
#' # the grid is adapted to the latest estimate, but weighted towards the prior 

init.quad <- function(Q = 2,
                      prior = list(mu = rep(0, Q), Sigma = diag(Q)),
                      adapt = NULL,
                      ip = 6,
                      prune = FALSE,
                      forcePD = FALSE,
                      debug = FALSE){
  
  # allow previous estimate input for adapt
  if (!is.null(adapt) && !is.null(attr(adapt, "variance"))){
    adapt <- list(mu = adapt, Sigma = attr(adapt, "variance"))
  }
  
  # make sure adapt is of the right size.
  if (!is.null(adapt) && (!is.list(adapt) || length(adapt$mu) != Q || dim(adapt$Sigma) != c(Q, Q))) stop("Size or format of Adapt argument invalid.")
  
  # get quadrature points, apply normal pdf
  x <- fastGHQuad::gaussHermiteData(ip)
  w <- x$w / sqrt(pi)
  x <- x$x * sqrt(2)
  
  # (if anyone knows an easy way to assign a single vector x times to x list elements, please tell me. )
  X <- as.matrix(expand.grid(lapply(apply(replicate(Q,x),2,list),unlist)))
  
  # compute lambda (eigen decomposition covar matrix)
  trans <- function(X, Sigma) {
    if(forcePD)
      Sigma <- nearPD(Sigma)$mat
    lambda <- with(eigen(Sigma), {
      if (any(values < 0)) 
        warning("Matrix is not positive definite.")
      if(length(values) > 1) 
        vectors %*% diag(sqrt(values))
      else 
        vectors * sqrt(values)
    })
    t(lambda %*% t(X))
  }
  
  # debugging, plot noise under the prior distribution and initial grid.
  if (debug){
    noise <- rmvnorm(10000, prior$mu, prior$Sigma)
    plot(noise,
         col = 'grey', pch = 20, cex = .2,
         xlim = range(noise[,1]) * 1.1,
         ylim = range(noise[,2]) * 1.1)
    points(X, pch = 20)
  }
  
  
  # calculate weights
  # same as above, roundabout way to get the combn weights for each combination of quad points
  g <- as.matrix(expand.grid(lapply(apply(replicate(Q,w),2,list),unlist)))
  # combined weight is the product of the individual weights, use sum of logs to reduce inaccuracy
  W <- apply(g,1,function(x) sum(log(x)))
  
  if (debug){
    points(X, pch = 20, col = 'green', cex = exp(W) / max(exp(W)) * 5)
    cat("\ngreen is transformed to prior")
  }
  
  # apply mv normal pdf error function if no earlier estimate is available
  if (is.null(adapt)){
    X <- trans(X, prior$Sigma)
    X <- t(t(X) + prior$mu)
  }
  # adapt to best estimate if there is
  else {
    # Adapt quadrature grid
    X <- trans(X, adapt$Sigma)
    X <- t(t(X) + adapt$mu)
    
    if (debug){
      points(X, pch = 20, col = 'blue', cex = exp(W) / max(exp(W)) * 5)
      cat("\nblue is transformed to posterior")
    }
    
    # calculate 'adaptive factor'
    # Ripped from mvtnorm/dmvnorm, note that the -.5 * Q * log(2*pi) term in each ll cancels out, the remainder is aux.
    adapt$chol <- chol(adapt$Sigma)
    adapt$det <- sum(log(diag(adapt$chol)))
    adapt$aux <- colSums(backsolve(adapt$chol, t(X) - adapt$mu, transpose = TRUE)^2)
    
    prior$chol <- chol(prior$Sigma)
    prior$det <- sum(log(diag(prior$chol)))
    prior$aux <- colSums(backsolve(prior$chol, t(X) - prior$mu, transpose = TRUE)^2)
    
    fact <- (adapt$aux - prior$aux) / 2 + adapt$det - prior$det
    
    # incorporate into weights.
    W <- W + fact
  }
  
  if(debug){
    points(X, pch = 20, col = 'orange', cex = exp(W) / max(exp(W)) * 5)
    cat('\norange is mean shifted')
  }
  
  # pruning
  if (prune) {
    threshold <- log(min(w) ^ (Q-1) * max(w))
    relevant <- W >= threshold
    
    W <- W[relevant]
    X <- X[relevant,,drop = FALSE]
    
    if(debug){
      points(X, pch = 20, col = 'pink', cex = exp(W) / max(exp(W)) * 5)
      cat('\npink is pruned.')
    }
  }
  
  
  return(invisible(list(X=X,W=W)))
}