#' Categorical Distribution
#'
#' Sample from a matrix of categorical probabilities.
#'
#' @param n Number of desired draws
#' @param p A \eqn{k}-dimensional vector of probabilities
#' @param log_p If `TRUE`, interpret given probabilities as having been
#' supplied on the log-scale. Otherwise, interpret them on the original scale.
#'
#' @returns a vector of category indices whose elements are in
#' \eqn{1, \ldots, k}.
#'
#' @details
#' We make use of the Gumbel trick to draw from probabilities given on the
#' log-scale without having to normalize. Note that `r_categ` can be slow
#' with large `n` because it runs in a loop in plain R.
#'
#' @name categ
NULL

#' @name categ
#' @export
r_categ = function(n, p, log_p = FALSE)
{
	k = length(p)
	lp = switch(log_p, "TRUE" = p, log(p))
	out = numeric(n)

	for (i in 1:n) {
		z = r_gumbel(k)
		out[i] = which.max(z + lp)
	}

	return(out)
}
