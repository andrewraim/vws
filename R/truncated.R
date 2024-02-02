#' Univariate Distributions Truncated to an Interval
#'
#' Functions to truncate a given univariate distribution to an interval of the
#' form \eqn{(\text{lo}, \text{hi}]}.
#'
#' @param n Desired sample size.
#' @param x Quantile or argument of density.
#' @param log logical; if \code{TRUE}, probabilities \code{p} are given as
#' \eqn{\log(p)}.
#' @param log.p logical; if \code{TRUE}, probabilities \code{p} are given as
#' \eqn{\log(p)}.
#' @param p Probability.
#' @param lo Lower limit.
#' @param hi Upper limit.
#' @param df Density function for the untruncated distribution.
#' @param pf Cumulative distribution function for the untruncated distribution.
#' @param qf Quantile function for the untruncated distribution.
#' @param ... Additional arguments.
#'
#' @details
#' This code was initially copied from the \link[LearnBayes]{rtruncated}
#' function in the LearnBayes package.
#'
#' @examples
#' # Take a large sample from truncated standard normal
#' x = r_truncated(n = 100000, lo = 1, hi = 2, pf = pnorm, qf = qnorm)
#'
#' # Compare density to histogram of draws
#' hist(x, probability = TRUE)
#' curve(d_truncated(x, lo = 1, hi = 2, df = dnorm, pf = pnorm), add = TRUE)
#'
#' # Compare CDF to empirical CDF of draws
#' plot(ecdf(x))
#' curve(p_truncated(x, lo = 1, hi = 2, pf = pnorm), col = "blue", lwd = 2, add = TRUE)
#'
#' # Compare quantile function to quantiles of draws
#' pr = seq(0, 1, length.out = 50)
#' q_emp = quantile(x, probs = pr)
#' q_thry = q_truncated(pr, lo = 1, hi = 2, pnorm, qf = qnorm)
#' plot(q_thry, q_emp)
#' abline(c(0,1), lty = 2, col = "red")
#'
#' @name TruncatedUnivariate
NULL

#' @name TruncatedUnivariate
#' @export
d_truncated = function(x, lo, hi, df, pf, log = FALSE, ...)
{
	log_p_lo = pf(lo, log.p = TRUE, ...)
	log_p_hi = pf(hi, log.p = TRUE, ...)
	out = df(x, log = TRUE, ...) + log(lo < x & x <= hi) - log_sub2_exp(log_p_hi, log_p_lo)
	if (log) { return(out) } else { return(exp(out)) }
}

#' @name TruncatedUnivariate
#' @export
p_truncated = function(x, lo, hi, pf, log.p = FALSE, ...)
{
	n = length(x)
	log_p_lo = pf(lo, log.p = TRUE, ...)
	log_p_hi = pf(hi, log.p = TRUE, ...)
	idx0 = which(x <= lo)
	idx1 = which(x > hi)
	idx2 = setdiff(1:n, c(idx0, idx1))
	out = numeric(n)
	out[idx0] = -Inf
	out[idx1] = 0
	log_p_x = pf(x[idx2], log.p = TRUE, ...)
	out[idx2] = log_sub2_exp(log_p_x, log_p_lo) - log_sub2_exp(log_p_hi, log_p_lo)
	if (log.p) { return(out) } else { return(exp(out)) }
}

#' @name TruncatedUnivariate
#' @export
q_truncated = function(p, lo, hi, pf, qf, log.p = FALSE, ...)
{
	if (log.p) { lp = p } else { lp = log(p) }
	log_p_lo = pf(lo, log.p = TRUE, ...)
	log_p_hi = pf(hi, log.p = TRUE, ...)
	log_p_adj = log_add2_exp(log_p_lo, lp + log_sub2_exp(log_p_hi, log_p_lo))
	qf(log_p_adj, log.p = TRUE, ...)
}

#' @name TruncatedUnivariate
#' @export
r_truncated = function(n, lo, hi, pf, qf, ...)
{
	u = runif(n)
	q_truncated(log(u), lo, hi, pf, qf, log.p = TRUE, ...)
}

