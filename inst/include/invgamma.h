#ifndef INVGAMMA_H
#define INVGAMMA_H

#include <Rcpp.h>
#include "log-sum-exp.h"

namespace vws {

//' Inverse Gamma distribution
//'
//' @param n Number of observations.
//' @param x Vector of quantiles.
//' @param q Vector of quantiles.
//' @param p Vector of probabilities.
//' @param a Shape parameter.
//' @param b Rate parameter.
//' @param lower logical; if TRUE (default), probabilities are
//' \eqn{P(X \leq x)}; otherwise, \eqn{P(X > x)}.
//' @param log If \code{TRUE}, return densities and probabilities on the log-scale.
//'
//' @return
//' \code{dinvgamma} gives the density, \code{rinvgamma} generates random
//' deviates.
//' @name InverseGamma
Rcpp::NumericVector r_invgamma(unsigned int n, double a, double b)
{
	return 1 / Rcpp::rgamma(n, a, 1/b);
}

//' @name InverseGamma
//' @export
double d_invgamma(double x, double a, double b, bool log = false)
{
	double out = R::dgamma(1/x, a, 1/b, true) - 2 * std::log(x);
	return log ? out : exp(out);
}

//' @name InverseGamma
//' @export
double p_invgamma(double q, double a, double b, bool lower = true, bool log = false)
{
	return R::pgamma(1 / q, a, 1/b, !lower, log);
}

//' @name InverseGamma
//' @export
double q_invgamma(double p, double a, double b, bool lower = true, bool log = false)
{
	if (!log) { p = std::log(p); }

	double out;
	if (p > std::log(1/2)) {
		double cp = log_sub2_exp(0, p);
		out = 1 / R::qgamma(cp, a, 1/b, lower, true);
	} else {
		out = 1 / R::qgamma(p, a, 1/b, !lower, true);
	}

	return out;
}

}

#endif
