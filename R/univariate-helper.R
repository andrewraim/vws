#' Univariate Distribution Helper
#'
#' A container to support VWS with univariate distributions.
#'
#' @param d Density function.
#' @param p Cumulative distribution function.
#' @param q Quantile function.
#'
#' @details
#' The specified functions must support the following interfaces.
#' \describe{
#' \item{`d(x, log = FALSE)`}{}
#' \item{`p(q, lower.tail = TRUE, log.p = FALSE)`}{}
#' \item{`q(p, lower.tail = TRUE, log.p = FALSE)`}{}
#' }
#'
#' Arguments to these functions are interpreted as usual for the `stats`
#' package. (Other arguments will be ignored if present).
#' \describe{
#' \item{`x`, `q`}{vector of quantiles.}{}
#' \item{`log`, `log.p`}{logical; if `TRUE`, probabilities
#' `p` are given as \eqn{\log(p)}.}{}
#' \item{`lower.tail`}{logical; if `TRUE`, probabilities are
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

	out = list(d = d, p = p, q = q)
	structure(out, class = "univariate_helper")
}

