d_bessel = function(x, a, nu, normalize = T, log = F)
{
	lnc = ifelse(normalize, log(besselI(a, nu)), 0)
	out = (2*x + nu) * log(a/2) - lfactorial(x) - lgamma(x + nu + 1) - lnc
	if (log) { return(out) } else { return(exp(out)) }
}
