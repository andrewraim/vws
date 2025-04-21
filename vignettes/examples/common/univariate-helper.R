#' Normal Distribution Helper
#'
#' A distribution helper based on \eqn{\text{N}(\mu, \sigma^2)}.
#'
#' @param mean Mean parameter \eqn{\mu}.
#' @param sd Standard deviation parameter \eqn{\sigma}.
#'
#' @examples
#' helper = normal_helper(mean = 0, sd = 1)
#' helper$d(0)
#' helper$p(1.96)
#' helper$q(0.025)
#'
#' @name normal_helper
#' @export
normal_helper = function(mean, sd)
{
	univariate_helper(
		d = function(x, log = FALSE) {
			dnorm(x, mean = mean, sd = sd, log = log)
		},
		p = function(q, lower.tail = TRUE, log.p = FALSE) {
			pnorm(q, mean = mean, sd = sd, lower.tail = lower.tail, log.p = log.p)
		},
		q = function(p, lower.tail = TRUE, log.p = FALSE) {
			qnorm(p, mean = mean, sd = sd, lower.tail = lower.tail, log.p = log.p)
		},
		s = function(x) { is.numeric(x) }
	)
}

#' Poisson Distribution Helper
#'
#' A distribution helper based on \eqn{\text{Poisson}(\lambda)}.
#'
#' @param lambda Rate parameter \eqn{\mu}.
#'
#' @examples
#' helper = poisson_helper(lambda = 10)
#' helper$d(5)
#' helper$p(10)
#' helper$q(0.025)
#'
#' @name poisson_helper
#' @export
poisson_helper = function(lambda)
{
	univariate_helper(
		d = function(x, log = FALSE) {
			dpois(x, lambda = lambda, log = log)}
		,
		p = function(q, lower.tail = TRUE, log.p = FALSE) {
			ppois(q, lambda = lambda, lower.tail = lower.tail, log.p = log.p)
		},
		q = function(p, lower.tail = TRUE, log.p = FALSE) {
			qpois(p, lambda = lambda, lower.tail = lower.tail, log.p = log.p)
		},
		s = function(x) { is.integer(x) & x >= 0 }
	)
}

#' Uniform Distribution Helper
#'
#' A distribution helper based on \eqn{\text{Uniform}(a, b)}.
#'
#' @param a lower limit of support.
#' @param b upper limit of support.
#'
#' @examples
#' helper = uniform_helper(a = 0, b = 1)
#' helper$d(1/2)
#' helper$p(1/4)
#' helper$q(0.025)
#'
#' @name uniform_helper
#' @export
uniform_helper = function(a, b)
{
	univariate_helper(
		d = function(x, log = FALSE) {
			dunif(x, min = a, max = b, log = log)
		},
		p = function(q, lower.tail = TRUE, log.p = FALSE) {
			punif(q, min = a, max = b, lower.tail = lower.tail, log.p = log.p)
		},
		q = function(p, lower.tail = TRUE, log.p = FALSE) {
			qunif(p, min = a, max = b, lower.tail = lower.tail, log.p = log.p)
		},
		s = function(x) { x >= a & x <= b }
	)
}

#' Inverse Gamma Distribution Helper
#'
#' A distribution helper based on \eqn{\text{IG}(a,b)}.
#'
#' @param a Shape parameter.
#' @param b Rate parameter.
#'
#' @examples
#' helper = invgamma_helper(a = 10, b = 5)
#' helper$d(1/2)
#' helper$p(0.29266)
#' helper$q(0.025)
#'
#' @name invgamma_helper
#' @export
invgamma_helper = function(a, b)
{
	univariate_helper(
		d = function(x, log = FALSE) {
			d_invgamma(x, a = a, b = b, log = log)
		},
		p = function(q, lower.tail = TRUE, log.p = FALSE) {
			p_invgamma(q, a = a, b = b, lower = lower.tail, log = log.p)
		},
		q = function(p, lower.tail = TRUE, log.p = FALSE) {
			q_invgamma(p, a = a, b = b, lower = lower.tail, log = log.p)
		},
		s = function(x) { x >= 0 }
	)
}


#' Gamma Distribution Helper
#'
#' A distribution helper based on \eqn{\text{Gamma}(a,b)}.
#'
#' @param a Shape parameter.
#' @param b Rate parameter.
#'
#' @examples
#' helper = gamma_helper(a = 10, b = 5)
#' helper$d(1/2)
#' helper$p(0.29266)
#' helper$q(0.025)
#'
#' @name invgamma_helper
#' @export
gamma_helper = function(a, b)
{
	univariate_helper(
		d = function(x, log = FALSE) {
			dgamma(x, shape = a, rate = b, log = log)
		},
		p = function(q, lower.tail = TRUE, log.p = FALSE) {
			pgamma(q, shape = a, rate = b, lower.tail = lower.tail, log.p = log.p)
		},
		q = function(p, lower.tail = TRUE, log.p = FALSE) {
			qgamma(p, shape = a, rate = b, lower.tail = lower.tail, log.p = log.p)
		},
		s = function(x) { x >= 0 }
	)
}
