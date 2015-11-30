library(MultiGHQuad)
require(Utilities)
Require("mvtnorm")

context("Results validate with naive grid quadrature")

pos.def <- function(n, ev = runif(n, 0, 10)) 
{
  if (n == 1) return(matrix(runif(1)))
  Z <- matrix(ncol = n, rnorm(n^2))
  decomp <- qr(Z)
  Q <- qr.Q(decomp) 
  R <- qr.R(decomp)
  d <- diag(R)
  ph <- d / abs(d)
  O <- Q %*% diag(ph)
  Z <- t(O) %*% diag(ev) %*% O
  return(Z)
}

# transform quad points to prior.
trans <- function(X, Sigma) {
  lambda <- with(eigen(Sigma), {
    if (any(values < 0)) warning("Matrix is not positive definite.")
    if(length(values) > 1) vectors %*% diag(sqrt(values))
    else vectors * sqrt(values)
  })
  t(lambda %*% t(X))
}

naive.quad <- function(prior, nQ, FUNC, ...){
  ## generate a naive grid with nQ equidistant nodes per dimension, centered on mu.
  nodes <- seq(-1, 1, length.out = nQ)
  plot(rmvnorm(5000, prior$mu, prior$Sigma), ylim = c(-3,3), xlim = c(-3,3), pch = 20, cex = .5, col = 'grey')
  
  # 2D
  grid <- expand.grid(nodes, nodes)
  points(grid)
  
  # scale to prior
  grid <- as.matrix(grid) %*% (4 * sqrt(prior$Sigma))
  points(grid, col = 'orange')
  
  # center at mu (rotate to add mu correctly).
  grid <- t(t(grid) + prior$mu)
  points(grid, col = 'blue')
  
  # evaluate the function
  evals <- apply(grid, 1, function(x) FUNC(x, ...))
  
  # normalizing constant
  evals <- evals / sum(evals)
  print(evals)
  print(grid)
  
  points(grid, cex = evals/max(evals) * 2, pch=20)
  print(colSums(evals * grid))
  invisible(list(evals, grid))
}

prior <- list(mu = c(2,1), Sigma = diag(2)/2)

naive.quad(prior, 6, dmvnorm, mean = prior$mu, sigma = prior$Sigma)
quadpts <- init.quad(Q = 2, prior = prior, ip = 10, adapt = prior, debug = TRUE)
eval.quad(dmvnorm, quadpts)
