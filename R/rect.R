#' Rectangular transformation
#'
#' A transformation from Euclidean to a rectangle in \eqn{\mathbb{R}^n}, and
#' its inverse transformation.
#'
#' @param x A point in \eqn{\mathbb{R}^{d}}.
#' @param z A point in the rectangle \eqn{[a_1,b_1] \times \cdots \times [a_d,b_d]}.
#' @param a A vector \eqn{(a_1, \ldots, a_d)}, Elements may be \code{-Inf}.
#' @param b A vector \eqn{(b_1, \ldots, b_d)}, Elements may be \code{Inf}.
#'
#' @name bounded
#' @examples
#' x = seq(-1, 1, length.out = 3)
#' a = rep(0, 3)
#' b = rep(1, 3)
#' z = inv_rect(x, a, b)
#' rect(z, a, b)
#'
#' a = c(-Inf, 0, -Inf)
#' b = c(Inf, 1, Inf)
#' z = inv_rect(x, a, b)
#' rect(z, a, b)
NULL

#' @name bounded
# @export
rect = function(z, a, b)
{
	n = length(x)
	stopifnot(n == length(a))
	stopifnot(n == length(b))
	stopifnot(all(a <= b))

	x = numeric(n)

	idx = which(is.infinite(a) & is.infinite(b) & a < 0 & b > 0)
	x[idx] = z[idx]

	idx = which(is.infinite(a) & a < 0 & is.finite(b))
	x[idx] = qlogis(z[idx] / b[idx])

	idx = which(is.infinite(b) & b > 0 & is.finite(a))
	x[idx] = log(z[idx] - a[idx])

	idx = which(is.finite(a) & is.finite(b))
	x[idx] = qlogis(z[idx] / (b[idx] - a[idx]) - a[idx])

	return(x)
}

#' @name bounded
# @export
inv_rect = function(x, a, b)
{
	n = length(x)
	stopifnot(n == length(a))
	stopifnot(n == length(b))
	stopifnot(all(a <= b))

	z = numeric(n)

	idx = which(is.infinite(a) & is.infinite(b) & a < 0 & b > 0)
	z[idx] = x[idx]

	idx = which(is.infinite(a) & a < 0 & is.finite(b))
	z[idx] = b[idx] * plogis(x[idx])

	idx = which(is.infinite(b) & b > 0 & is.finite(a))
	z[idx] = exp(x[idx]) + a[idx]

	idx = which(is.finite(a) & is.finite(b))
	z[idx] = (b[idx] - a[idx]) * plogis(x[idx]) + a[idx]

	return(z)
}

