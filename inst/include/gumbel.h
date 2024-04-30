#ifndef GUMBEL_H
#define GUMBEL_H

#include <Rcpp.h>

namespace vws {

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
Rcpp::NumericVector q_gumbel(const Rcpp::NumericVector& p, double mu = 0,
	double sigma = 1, bool lower = true, bool log = false)
{
	// const Rcpp::NumericVector& lp0 = log ? p : Rcpp::log(p);
	Rcpp::NumericVector lp0;
	if (log) {
		lp0 = p;
	} else {
		lp0 = Rcpp::log(p);
	}

	const Rcpp::NumericVector& lp = lower ? lp0 : Rcpp::log1p(-Rcpp::exp(lp0));
	return mu - sigma * Rcpp::log(-lp);
}

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
	const Rcpp::NumericVector& z = (x - mu) / sigma;
	const Rcpp::NumericVector& out = -std::log(sigma) - (z + Rcpp::exp(-z));
	if (log) {
		return out;
	} else {
		return Rcpp::exp(out);
	}
}

//' @name Gumbel
//' @export
Rcpp::NumericVector p_gumbel(const Rcpp::NumericVector& q, double mu = 0,
	double sigma = 1, bool lower = true, bool log = false)
{
	const Rcpp::NumericVector& z = (q - mu) / sigma;
	const Rcpp::NumericVector& out0 = -exp(-z);
	const Rcpp::NumericVector& out = lower ? out0 : Rcpp::log1p(-Rcpp::exp(out0));
	if (log) {
		return out;
	} else {
		return Rcpp::exp(out);
	}
}

}

#endif
