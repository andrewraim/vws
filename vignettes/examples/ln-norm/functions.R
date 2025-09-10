library(vws)
library(statmod)

d_target = function(x, mu, sigma, z, lambda, log = FALSE, normalize = TRUE, Q = 10)
{
	psi = 1
	if (normalize) {
		## Compute normalizing constant using Hermite quadrature with Q points
		w = function(x, log = TRUE) {
			out = -log(x) - (log(x) - mu)^2 / (2*sigma^2) + log(x > 0)
			out[x == 0] = -Inf
			if (log) { return(out) } else { return(exp(out)) }
		}

		q = function(x, log = FALSE) {
			tx = z - sqrt(2) * lambda * x
			ii = which(tx > 0)
			out = rep(-Inf, length(x))
			out[ii] = -log(tx[ii]) - (log(tx[ii]) - mu)^2 / (2*sigma^2) - 1/2 * log(pi)
			if (log) { return(out) } else { return(exp(out)) }
		}

		quad_out = gauss.quad(n = Q, kind = "hermite")
		ww = quad_out$weights
		xx = quad_out$nodes
		psi = sum(ww * q(xx))
	}

	out = w(x, log = TRUE) + dnorm(x, z, lambda, log = TRUE) - log(psi)
	if (log) { return(out) } else { return(exp(out)) }
}
