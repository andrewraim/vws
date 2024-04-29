#include <Rcpp.h>

//' Gumbel Distribution
//'
//' Functions for the Gumbel distribution
//'
//' @param n Number of desired draws
//' @param x Vector of quantiles
//' @param p Vector of probabilities
//' @param q Vector of quantiles
//' @param mu Location parameter
//' @param sigma Scale parameter
//' @param lower Logical; if \code{TRUE} (default), probabilities are
//' \eqn{P[X \leq x]} otherwise, \eqn{P[X > x]}.
//' @param log Logical; if \code{TRUE}, probabilities p are given as \eqn{log(p)}
//'
//' @return A vector of draws
//'
//' @name Gumbel
NULL

//' @name Gumbel
//' @export
Rcpp::NumericVector r_gumbel(unsigned int n, double mu = 0, double sigma = 1)
{
	const Rcpp::NumericVector& u = Rcpp::runif(n);
	return q_gumbel(u, mu, sigma);
}

//' @name Gumbel
//' @export
Rcpp::NumericVector d_gumbel(const Rcpp::NumericVector& x, double mu = 0,
	double sigma = 1, bool log = false)
{
	z = (x - mu) / sigma;
	out = -log(sigma) - (z + exp(-z));
	return log ? out : exp(out);
}

//' @name Gumbel
//' @export
Rcpp::NumericVector p_gumbel(const Rcpp::NumericVector& q, double mu = 0,
	double sigma = 1, bool lower = true, log = false)
{
	z = (q - mu) / sigma;
	out0 = -exp(-z);
	if (lower) { out = out0; } else { out = log1p(-exp(out0)); }
	return log ? out : exp(out);
}

//' @name Gumbel
//' @export
Rcpp::NumericVector q_gumbel(const Rcpp::NumericVector& p, double mu = 0,
	double sigma = 1, bool lower = true, bool log = false)
{
	if (log) { lp0 = p; } else { lp0 = log(p); }
	if (lower) { lp = lp0; } else { lp = log1p(-exp(lp0)); }
	return mu - sigma * log(-lp);
}
