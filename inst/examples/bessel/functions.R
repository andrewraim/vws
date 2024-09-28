d_bessel = function(x, a, nu, normalize = TRUE, log = FALSE)
{
	if (normalize) {
		out = (2*x + nu) * log(a/2) - log(besselI(a, nu)) - lfactorial(x) - lgamma(x + nu + 1)
	} else {
		# out = (2*x) * log(a/2) - lfactorial(x) - lgamma(x + nu + 1)
		# out = x*log(a^2/4) - a^2/4 - lfactorial(x) - lgamma(x + nu + 1)
		out = (2*x + nu) * log(a/2) - lfactorial(x) - lgamma(x + nu + 1)
	}
	if (log) { return(out) } else { return(exp(out)) }
}
