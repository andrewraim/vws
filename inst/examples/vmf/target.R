# Operations for Target Distribution

## Normalizing constant.
nc_target = function(kappa, d, log = FALSE)
{
	out = 0.5*log(base::pi) - (d/2 - 1)*log(kappa / 2) +
		log(besselI(kappa, d/2 - 1)) + lgamma(d/2 - 0.5) +
		log(kappa) - vws::log_sub2_exp(kappa,-kappa)
	if (log) { return(out) } else { return(exp(out)) }
}

## Density.
d_target = function(x, kappa, d, normalize = TRUE, log = FALSE)
{
	n = length(x)
	out = rep(-Inf, n)
	idx = which(-1 < x & x < 1)
	out[idx] = 0.5 * (d-3) * log1p(-x[idx]^2) +
		d_base(x[idx], kappa, log = TRUE)
	if (normalize) { out = out - nc_target(kappa, d, log = TRUE) }
	if (log) { return(out) } else { return(exp(out)) }
}

## CDF.
p_target = function(x, kappa, d, lower.tail = TRUE, log.p = FALSE)
{
	n = length(x)
	islower = rep(lower.tail, n)
	lo = ifelse(islower, rep(-1, n), x)
	hi = ifelse(!islower, rep(1, n), x)

	out = numeric(n)
	for (i in 1:n) {
		if (lo[i] < hi[i]) {
			int_out = integrate(d_target, lower = lo[i], upper = hi[i],
				kappa = kappa, d = d)
			out[i] = int_out$value
		}
	}
	if (log.p) { return(log(out)) } else { return(out) }
}
