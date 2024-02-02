#' Gumbel Distribution
#' 
#' Functions for the Gumbel distribution
#' 
#' @param n Number of desired draws
#' @param x Vector of quantiles
#' @param p Vector of probabilities
#' @param q Vector of quantiles
#' @param mu Location parameter
#' @param sigma Scale parameter
#' @param lower.tail Logical; if \code{TRUE} (default), probabilities are
#' \eqn{P[X \leq x]} otherwise, \eqn{P[X > x]}.
#' @param log.p Logical; if \code{TRUE}, probabilities p are given as \eqn{log(p)}
#' @param log Logical; if \code{TRUE}, probabilities p are given as \eqn{log(p)}
#' 
#' @return A vector of draws
#' 
#' @name Gumbel
NULL

#' @name Gumbel
#' @export
r_gumbel = function(n, mu = 0, sigma = 1)
{
	u = runif(n)
	q_gumbel(u, mu, sigma)
}

#' @name Gumbel
#' @export
d_gumbel = function(x, mu = 0, sigma = 1, log = FALSE)
{
	z = (x - mu) / sigma
	out = -log(sigma) - (z + exp(-z))
	if (log) { return(out) } else { return(exp(out)) }
}

#' @name Gumbel
#' @export
p_gumbel = function(q, mu = 0, sigma = 1, lower.tail = TRUE, log.p = FALSE) 
{
	z = (q - mu) / sigma
	out0 = -exp(-z)
	if (lower.tail) { out = out0 } else { out = log1p(-exp(out0)) }
	if (log.p) { return(out) } else { return(exp(out)) }
}

#' @name Gumbel
#' @export
q_gumbel = function(p, mu = 0, sigma = 1, lower.tail = TRUE, log.p = FALSE) 
{
	if (log.p) { lp0 = p } else { lp0 = log(p) }
	if (lower.tail) { lp = lp0 } else { lp = log1p(-exp(lp0)) }
	mu - sigma * log(-lp)
}

