r_mvnorm_prec = function(mu, Omega)
{
	stopifnot(length(mu) == nrow(Omega))
	stopifnot(nrow(Omega) == ncol(Omega))

    k = length(mu)
    z = rnorm(k)
    A = chol(Omega)
    out = solve(A, z) + mu

	return(out)
}
