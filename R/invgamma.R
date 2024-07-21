#' Inverse Gamma distribution
#'
#' @param n Number of observations.
#' @param x Vector of quantiles.
#' @param q Vector of quantiles.
#' @param p Vector of probabilities.
#' @param a Shape parameter.
#' @param b Rate parameter.
#' @param lower.tail logical; if TRUE (default), probabilities are
#' \eqn{P(X \leq x)}; otherwise, \eqn{P(X > x)}.
#' @param log If \code{TRUE}, return densities and probabilities on the log-scale.
#' @param log.p If \code{TRUE}, input and output probabilities on the log-scale.
#'
#' @return
#' \code{dinvgamma} gives the density, \code{rinvgamma} generates random
#' deviates.
#' @name InverseGamma
NULL

#' @name InverseGamma
#' @export
r_invgamma = function(n, a, b)
{
	1 / rgamma(n, a, b)
}

#' @name InverseGamma
#' @export
d_invgamma = function(x, a, b, log = FALSE)
{
	out = dgamma(1/x, a, b, log = TRUE) - 2 * log(x)
	if (log) { return(out) } else { return(exp(out))}
}

#' @name InverseGamma
#' @export
p_invgamma = function(q, a, b, lower.tail = TRUE, log.p = FALSE)
{
	pgamma(1 / q, a, b, lower.tail = !lower.tail, log.p = log.p)
}

#' @name InverseGamma
#' @export
q_invgamma = function(p, a, b, lower.tail = TRUE, log.p = FALSE)
{
	n = length(p)
	if (!log.p) { p = log(p) }
	if (length(a) == 1) { a = rep(a, n)}
	if (length(b) == 1) { b = rep(b, n)}

	out = numeric(n)
	cp = log_sub2_exp(0, p)

	for (i in 1:n) {
		if (p[i] > log(1/2)) {
			out[i] = 1 / qgamma(cp[i], a[i], b[i], lower.tail = lower.tail, log.p = TRUE)
		} else {
			out[i] = 1 / qgamma(p[i], a[i], b[i], lower.tail = !lower.tail, log.p = TRUE)
		}
	}

	return(out)
}

