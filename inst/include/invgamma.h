#ifndef VWS_INVGAMMA_H
#define VWS_INVGAMMA_H

#include <Rcpp.h>
#include "log-sum-exp.h"

namespace vws {

/*
* Non-vectorized distribution functions.
*/

inline double r_invgamma(double a, double b)
{
	return 1 / R::rgamma(a, 1 / b);
}

inline double d_invgamma(double x, double a, double b, bool log = false)
{
	double out = R::dgamma(1 / x, a, 1 / b, true) - 2 * std::log(x);
	return log ? out : exp(out);
}

inline double p_invgamma(double q, double a, double b, bool lower = true, bool log = false)
{
	return R::pgamma(1 / q, a, 1 / b, !lower, log);
}

inline double q_invgamma(double p, double a, double b, bool lower = true, bool log = false)
{
	if (!log) { p = std::log(p); }

	double out;
	if (p > std::log(1/2)) {
		double cp = log_sub2_exp(0, p);
		out = 1 / R::qgamma(cp, a, 1 / b, lower, true);
	} else {
		out = 1 / R::qgamma(p, a, 1 / b, !lower, true);
	}

	return out;
}

/*
* Vectorized distribution functions.
*/

inline Rcpp::NumericVector r_invgamma(unsigned int n, double a, double b)
{
	return 1 / Rcpp::rgamma(n, a, 1 / b);
}

inline Rcpp::NumericVector d_invgamma(const Rcpp::NumericVector& x, double a,
	double b, bool log)
{
	unsigned int n = x.size();
	Rcpp::NumericVector out(n);

	for (unsigned int i = 0; i < n; i++) {
		out(i) = vws::d_invgamma(x(i), a, b, log);
	}

	return out;
}

inline Rcpp::NumericVector p_invgamma(const Rcpp::NumericVector& q, double a,
	double b, bool lower, bool log)
{
	unsigned int n = q.size();
	Rcpp::NumericVector out(n);

	for (unsigned int i = 0; i < n; i++) {
		out(i) = vws::p_invgamma(q(i), a, b, log);
	}

	return out;
}

inline Rcpp::NumericVector q_invgamma(const Rcpp::NumericVector& p, double a,
	double b, bool lower, bool log)
{
	unsigned int n = p.size();
	Rcpp::NumericVector out(n);

	for (unsigned int i = 0; i < n; i++) {
		out(i) = vws::q_invgamma(p(i), a, b, log);
	}

	return out;
}

}

#endif
