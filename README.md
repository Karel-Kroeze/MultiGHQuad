MultiGHQuad
===========
*Multidimensional Gauss-Hermite quadrature in R.*

The package consists of two functions, eval.quad() and init.quad(). The former evaluates the integral of any given function with respect to Q multivariate normal distributed parameters. The latter provides the quadrature points and weights at which the function should be evaluated, and the weights attached to each evaluation.

Quadrature points are taken from package fastGHQuad, and expanded into a multidimensional grid. This grid is then adjusted for a given multivariate normal distribution - usually a prior or population distribution. When this prior distribution is far from the posterior distribution under integration, precision of the estimate suffers. To overcome this, the quadrature grid can be further adjusted to adapt to an estimate of the posterior distribution.

&copy; Karel Kroeze 2014.

[![Build Status](https://travis-ci.org/Karel-Kroeze/MultiGHQuad.svg?branch=master)](https://travis-ci.org/Karel-Kroeze/MultiGHQuad)
