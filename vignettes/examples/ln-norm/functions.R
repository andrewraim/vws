library(vws)
library(statmod)

get_target_density = function(mu, sigma, z, lambda)
{
	w = function(y, log = TRUE) {
		out = -log(y) - (log(y) - mu)^2 / (2*sigma^2) + log(y > 0)
		out[y == 0] = -Inf
		if (log) { return(out) } else { return(exp(out)) }
	}
	q = function(x, log = FALSE) {
		tx = z - sqrt(2) * lambda * x
		ii = which(tx > 0)
		out = rep(-Inf, length(x))
		out[ii] = -log(tx[ii]) - (log(tx[ii]) - mu)^2 / (2*sigma^2) - 1/2 * log(pi)
		if (log) { return(out) } else { return(exp(out)) }
	}
	quad_out = gauss.quad(n = 10, kind = "hermite")
	ww = quad_out$weights
	xx = quad_out$nodes
	psi = sum(ww * q(xx))
	d_target = function(y, log = TRUE) {
		out = w(y, log = TRUE) + dnorm(y, z, lambda, log = TRUE) - log(psi)
		if (log) { return(out) } else { return(exp(out)) }
	}

	return(d_target)
}
