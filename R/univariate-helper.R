#' Univariate Distribution Helper
#'
#' A container to support VWS with univariate distributions.
#'
#' @param d Density function.
#' @param p Cumulative distribution function.
#' @param q Quantile function.
#' @param s Indicator function that returns \code{TRUE} when the
#' argument is in the support of the distribution; otherwise returns
#' \code{FALSE}.
#'
#' @details
#' The specified functions must support the following interfaces.
#' \describe{
#' \item{\code{d(x, log = FALSE)}}{}
#' \item{\code{p(q, lower.tail = TRUE, log.p = FALSE)}}{}
#' \item{\code{q(p, lower.tail = TRUE, log.p = FALSE)}}{}
#' \item{\code{s(x)}}{}
#' }
#'
#' Arguments to these functions are interpreted as usual for the \code{stats}
#' package. (Other arguments will be ignored if present).
#' \describe{
#' \item{\code{x}, \code{q}}{vector of quantiles.}{}
#' \item{\code{log}, \code{log.p}}{logical; if \code{TRUE}, probabilities
#' \code{p} are given as \eqn{\log(p)}.}{}
#' \item{\code{lower.tail}}{logical; if \code{TRUE}, probabilities are
#' \eqn{\text{P}(X \leq x)}; otherwise \eqn{\text{P}(X > x)}}{}
#' }
#'
#' @examples
#' helper = univariate_helper(
#'     d = function(x, log = FALSE) dpois(x, 10, log),
#'     p = function(q, lower.tail = TRUE, log.p = FALSE) {
#'         ppois(q, 10, lower.tail, log.p)
#'     },
#'     q = function(p, lower.tail = TRUE, log.p = FALSE) {
#'         qpois(p, 10, lower.tail, log.p)
#'     },
#'     s = function(x) { is.numeric(x) }
#' )
#'
#' helper$d(5)
#' helper$p(10)
#' helper$q(0.025)
#'
#' @name univariate_helper
#' @export
univariate_helper = function(d, p, q, s)
{
	# Make sure functions have the needed arguments
	stopifnot(all(c("x", "log") %in% names(formals(d))))
	stopifnot(all(c("q", "lower.tail", "log.p") %in% names(formals(p))))
	stopifnot(all(c("p", "lower.tail", "log.p") %in% names(formals(q))))
	stopifnot(all(c("x") %in% names(formals(s))))

	out = list(d = d, p = p, q = q, s = s)
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
