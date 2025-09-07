d_bessel = function(x, a, nu, normalize = TRUE, log = FALSE)
{
	lnc = ifelse(normalize, log(besselI(a, nu)), 0)
	out = (2*x + nu) * log(a/2) - lfactorial(x) - lgamma(x + nu + 1) - nc
	if (log) { return(out) } else { return(exp(out)) }
}

