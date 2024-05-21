library(statmod)

nc_target = function(mu, sigma2, z, lambda2, log = FALSE)
{
	q = function(x, log = FALSE) {
		tx = z - sqrt(2 * lambda2) * x
		out = log(tx > 0) - log(tx) - (log(tx) - mu)^2 / (2*sigma2) - 1/2 * log(pi)
		if (log) { return(out) } else { return(exp(out)) }
	}

	quad_out = gauss.quad(n = 10, kind = "hermite")
	ww = quad_out$weights
	xx = quad_out$nodes
	psi = sum(ww * q(xx))

	if (log) { return(log(psi)) } else { return(psi) }
}

d_target = function(y, mu, sigma2, z, lambda2, log = TRUE)
{
	log_wx = -log(y) - (log(y) - mu)^2 / (2*sigma2) + log(y > 0)
	log_wx[y == 0] = -Inf

	log_gx = dnorm(y, mean = z, sd = sqrt(lambda2), log = TRUE)

	log_nc = nc_target(mu, sigma2, z, lambda2, log = TRUE)

	out = log_wx + log_gx - log_nc
	if (log) { return(out) } else { return(exp(out)) }
}
