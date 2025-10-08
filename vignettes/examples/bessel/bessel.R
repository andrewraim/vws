w = function(x, log = TRUE) {
	out = -lgamma(x + nu + 1)
	if (log) { return(out) } else { return(exp(out)) }
}

d_bessel = function(x, a, nu, normalize = T, log = F)
{
	lnc = ifelse(normalize, log(besselI(a, nu)), 0)
	out = (2*x + nu) * log(a/2) - lfactorial(x) - lgamma(x + nu + 1) - lnc
	if (log) { return(out) } else { return(exp(out)) }
}

d_truncpois = function(x, lambda, a, b, log = F)
{
	lo = max(a, 0)
	hi = max(b, 0)
	# lg1 = incgamma(floor(hi) + 1, lambda, lower = F, log = T)
	# lg2 = incgamma(floor(hi) + 1, 0, lower = F, log = T)
	# lg3 = incgamma(ceiling(lo), lambda, lower = F, log = T)
	# lg4 = incgamma(ceiling(lo), 0, lower = F, log = T)
	# lp = vws::log_sub2_exp(lg1 - lg2, lg3 - lg4)

	lpa = ppois(lo, lambda, lower.tail = T, log.p = T)
	lpb = ppois(hi, lambda, lower.tail = T, log.p = T)
	clpa = ppois(lo, lambda, lower.tail = F, log.p = T)
	clpb = ppois(hi, lambda, lower.tail = F, log.p = T)
	lp = vws::log_sub2_exp(lpb, lpa)
	clp = vws::log_sub2_exp(clpa, clpb)
	lm = max(lp, clp)

	out = x * log(lambda) - lambda - lfactorial(x) - lm + log(lo < x & x <= hi)
	if (log) { return(out) } else { return(exp(out)) }
}
