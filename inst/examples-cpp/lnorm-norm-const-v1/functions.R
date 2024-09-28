library(vws)
library(statmod)

get_target_density = function(mu, sigma2, z, lambda2)
{
	w = function(y, log = TRUE) {
		out = -log(y) - (log(y) - mu)^2 / (2*sigma2) + log(y > 0)
		out[y == 0] = -Inf
		if (log) { return(out) } else { return(exp(out)) }
	}
	helper = normal_helper(mean = z, sd = sqrt(lambda2))
	q = function(x, log = FALSE) {
		tx = z - sqrt(2 * lambda2) * x
		out = log(tx > 0) - log(tx) - (log(tx) - mu)^2 / (2*sigma2) - 1/2 * log(pi)
		if (log) { return(out) } else { return(exp(out)) }
	}
	quad_out = gauss.quad(n = 10, kind = "hermite")
	ww = quad_out$weights
	xx = quad_out$nodes
	psi = sum(ww * q(xx))
	d_target = function(y, log = TRUE) {
		out = w(y, log = TRUE) + helper$d(y, log = TRUE) - log(psi)
		if (log) { return(out) } else { return(exp(out)) }
	}

	return(d_target)
}
