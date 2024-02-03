#' Univariate Distribution Helper
#'
#' A container to support VWS with univariate distributions.
#'
#' @param d Density function for base distribution.
#' @param p Cumulative distribution function for base distribution.
#' @param q Quantile function for base distribution.
#' @param in_support Support inclusion function for base distribution.
#'
#' @examples
#' univariate_helper(
#'     d = function(x, log = FALSE) dpois(x, log),
#'     p = function(q, lower.tail = TRUE, log.p = FALSE) {
#'         ppois(q, lower.tail, log.p)
#'     },
#'     q = function(p, lower.tail = TRUE, log.p = FALSE) {
#'         qpois(p, lower.tail, log.p)
#'     },
#'     in_support = function(x) { TRUE }
#' )
#'
#' @name univariate_helper
#' @export
univariate_helper = function(d, p, q, in_support)
{
	# Make sure functions have the needed arguments
	stopifnot(all(c("x", "log") %in% names(formals(d))))
	stopifnot(all(c("q", "lower.tail", "log.p") %in% names(formals(p))))
	stopifnot(all(c("p", "lower.tail", "log.p") %in% names(formals(q))))
	stopifnot(all(c("x") %in% names(formals(in_support))))

	out = list(d = d, p = p, q = q, in_support = in_support)
	structure(out, class = "univariate_helper")
}

#' Normal Distribution Helper
#'
#' A distribution helper based on \eqn{\text{N}(\mu, \sigma^2)}.
#'
#' @param mean Mean parameter \eqn{\mu}.
#' @param sd Standard deviation parameter \eqn{\sigma}.
#'
#' @examples
#' helper = normal_univariate_helper(mean = 0, sd = 1)
#' helper$d(0)
#' helper$p(1.96)
#' helper$q(0.025)
#'
#' @name normal_univariate_helper
#' @export
normal_univariate_helper = function(mean, sd)
{
	univariate_helper(
		d = function(x, log = FALSE) dnorm(x, mean, sd, log),
		p = function(q, lower.tail = TRUE, log.p = FALSE) {
			pnorm(q, mean, sd, lower.tail, log.p)
		},
		q = function(p, lower.tail = TRUE, log.p = FALSE) {
			qnorm(p, mean, sd, lower.tail, log.p)
		},
		in_support = function(x) { TRUE }
	)
}

#' Poisson Distribution Helper
#'
#' A distribution helper based on \eqn{\text{Poisson}(\lambda)}.
#'
#' @param lambda Rate parameter \eqn{\mu}.
#'
#' @examples
#' helper = poisson_univariate_helper(lambda = 10)
#' helper$d(5)
#' helper$p(10)
#' helper$q(0.025)
#'
#' @name poisson_univariate_helper
#' @export
poisson_univariate_helper = function(lambda)
{
	univariate_helper(
		d = function(x, log = FALSE) dpois(x, lambda, log),
		p = function(q, lower.tail = TRUE, log.p = FALSE) {
			ppois(q, lambda, lower.tail, log.p)
		},
		q = function(p, lower.tail = TRUE, log.p = FALSE) {
			qpois(p, lambda, lower.tail, log.p)
		},
		in_support = function(x) { is.integer(x) & x >= 0 }
	)
}
