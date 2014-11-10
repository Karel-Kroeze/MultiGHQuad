MultiGHQuad
===========
Multidimensional Gauss-Hermite quadrature in R.

The package consists of two functions, MGHQuadEval() and MGHQuadPoints(). The former evaluates the integral of any given function with respect to Q multivariate normal distributed parameters. The latter provides the quadrature points and weights at which the function should be evaluated, and the weights attached to each evaluation.

(c) Karel Kroeze 2014.

Quadrature points are taken from package fastGHQuad, and expanded into a multidimensional grid.
