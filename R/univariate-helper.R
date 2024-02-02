#' Univariate Distribution Helper
#'
#' TBD: document expected function interface.
#'
#' @param r Drawing function for base distribution.
#' @param d Density function for base distribution.
#' @param p CDF function for base distribution.
#' @param q Quantile function for base distribution.
#' @param in_support Support inclusion function for base distribution.
#'
#' @name univariate_helper
#' @export
univariate_helper = function(r, d, p, q, in_support)
{
	# Make sure functions have the needed arguments
	stopifnot(all(c("n") %in% names(formals(r))))
	stopifnot(all(c("x", "log") %in% names(formals(d))))
	stopifnot(all(c("q", "lower.tail", "log.p") %in% names(formals(p))))
	stopifnot(all(c("p", "lower.tail", "log.p") %in% names(formals(q))))
	stopifnot(all(c("x") %in% names(formals(in_support))))

	out = list(r = r, d = d, p = p, q = q, in_support = in_support)
	structure(out, class = "univariate_helper")
}
