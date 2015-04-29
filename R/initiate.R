#' Q-dimensional grid of quadrature points.
#' 
#' Creates a flattened, rotated grid that incorporates correlation through an eigenvalue decomposition of the covariance matrix.
#'
#' Creates a Q-dimensional grid by calling \code{\link{expand.grid}} on Q vectors of unidimensional quadrature points.
#' The grid is corrected for covariation by eigenvalue decomposition;
#' \deqn{\Sigma = S \times \Lambda \times S^T}{ Sigma = S \%*\% Lambda \%*\% t(S) }
#' Take \code{A} to be;
#' \deqn{A = S \times \sqrt(\Lambda)}{A = S \%*\% \sqrt(\Lambda)}
#' And left multiply the quadrature points \code{z} by \code{A} to obtain correlated quadrature points \code{r};
#' \deqn{r = A \times z}{r = A \%*\% z}
#' 
#' @param Q Number of dimensions. Defaults to 2. Only required when \code{mu} and \code{Sigma} are not provided.
#' @param prior List of prior mean \code{mu}, = \code{vector}, and covariance matrix \code{Sigma} = \code{matrix}, defaults to zero vector and identity matrix respectively.
#' @param adapt List of adaptive mean \code{mu}, = \code{vector}, and covariance matrix \code{Sigma} = \code{matrix}, if \code{NULL} no adaptation is used. Defaults to NULL.
#' @param ip Number of quadrature points \emph{per dimension}. Defaults to 6. Note that the total number of quadrature points is \code{ip^Q}.
#' @return A list with a matrix \code{X} of \code{ip^Q} by \code{Q} quadrature points and a vector \code{W} of length \code{ip^Q} associated weights.
#' @seealso \code{\link[fastGHQuad]{gaussHermiteData}}, used to create unidimensional quadrature points, and \code{\link{eval.quad}} for evaluating the integral.
#' @export
#' @examples
#' # generate some noise with a given covariance matrix
#' \dontrun{
#' require(mvtnorm)
#' sigma <- matrix(c(1,.8,.8,1),ncol=2,byrow=T)
#' noise <- rmvnorm(1e4,mean=c(0,0),sigma=sigma)
#' # plot noise
#' plot(noise,col='red',pch='.')
#' 
#' # generate quadrature points
#' quadPoints <- init.quad(prior = list(mu=c(0,0),Sigma=sigma),ip=10)
#' 
#' # plot quad points
#' points(quadPoints$X,pch=16)
#' 
#' # plot quad points with (log) weights
#' plot(noise,col='red',pch='.')
#' points(quadPoints$X,col=grey(1-quadPoints$W/max(quadPoints$W)),pch=16)
#' }
init.quad <- function(Q = 2,
                      prior = list(mu = rep(0, Q), Sigma = diag(Q)),
                      adapt = NULL,
                      ip=6){
  # TODO: account for mu != 0.
  # get quadrature points, apply normal pdf
  x <- fastGHQuad::gaussHermiteData(ip)
  w <- x$w / sqrt(pi)
  x <- x$x * sqrt(2)
  
  # (if anyone knows an easy way to assign a single vector x times to x list elements, please tell me. )
  X <- as.matrix(expand.grid(lapply(apply(replicate(Q,x),2,list),unlist)))
  
  # compute lambda (eigen decomposition covar matrix)
  prior$lambda <- with(eigen(prior$Sigma), vectors %*% diag(sqrt(values)))
  
  # apply mv normal pdf error function 
  # account for correlation
  X <- t(prior$lambda %*% t(X))
  # plug in prior mean
  # TODO: Check with Cees if it's really this simple
  X <- t(t(X) + prior$mu)
  
  # calculate weights
  # same as above, roundabout way to get the combn weights for each combination of quad points
  g <- as.matrix(expand.grid(lapply(apply(replicate(Q,w),2,list),unlist)))
  # combined weight is the product of the individual weights
  W <- log(apply(g,1,prod))
  
  # adapt to best estimate
  if (!is.null(adapt)){
    # Adapt quadrature grid
    adapt$lambda <- with(eigen(adapt$Sigma), vectors %*% diag(sqrt(values)))
    X <- t(adapt$lambda %*% t(X))
    X <- t(t(X) + adapt$mu)
    
    # calculate 'adaptive factor'
    # Ripped from mvtnorm/dmvnorm, note that the -.5 * Q * log(2*pi) term in each ll cancels out, the remainder is aux.
    adapt$chol <- chol(adapt$Sigma)
    adapt$det <- sum(log(diag(adapt$chol)))
    adapt$aux <- colSums(backsolve(adapt$chol, t(X) - adapt$mu, transpose = TRUE)^2)
    
    prior$chol <- chol(prior$Sigma)
    prior$det <- sum(log(diag(prior$chol)))
    prior$aux <- colSums(backsolve(prior$chol, t(X) - prior$mu, transpose = TRUE)^2)
    
    fact <- adapt$det - prior$det + (adapt$aux - prior$aux) / 2

    # incorporate into weights.
    W <- W + fact
  }
  
  return(invisible(list(X=X,W=W)))
}