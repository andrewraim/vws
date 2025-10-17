#ifndef TARGET_H
#define TARGET_H

#include "vws.h"

// [[Rcpp::export(name = "n_target")]]
inline double n_target(double kappa, double d, bool log = false)
{
	double bes = R::bessel_i(kappa, d / 2 - 1, 1);
	double out = 0.5 * std::log(M_PI) - (d / 2 - 1) * std::log(kappa / 2) +
		std::log(bes) + lgamma(d / 2 - 0.5);
	return log ? out : std::exp(out);
}

inline double d_target(double x, double kappa, double d, bool log = false,
	bool normalize = true)
{
	double log_nc = normalize ? n_target(kappa, d, true) : 0;

	double out = R_NegInf;
	if (-1 < x && x < 1) {
		out = 0.5 * (d - 3) * std::log1p(-std::pow(x, 2)) + kappa*x - log_nc;
	}

	return log ? out : std::exp(out);
}

inline double p_target(double q, double kappa, double d, bool lower = true,
	bool log = false)
{
	double lo = lower ? -1 : q;
	double hi = lower ? q : 1;

	const fntl::dfd& f = [&](double x) -> double {
		return d_target(x, kappa, d);
	};
	const auto& int_out = fntl::integrate(f, lo, hi);
	double out = int_out.value;

	return log ? std::log(out) : out;
}

// Compute integrals, from a to b, under the normalized target
// [[Rcpp::export(name = "integrate_target")]]
inline double integrate_target(double a, double b, double kappa, double d,
	bool log = false)
{
	const fntl::dfd& f = [&](double x) -> double {
		return d_target(x, kappa, d);
	};
	const auto& int_out = fntl::integrate(f, a, b);
	double out = std::log(int_out.value);
	return log ? out : std::exp(out);
}

// [[Rcpp::export]]
inline Rcpp::NumericVector d_target(const Rcpp::NumericVector& x, double kappa,
	double d, bool log = false, bool normalize = true)
{
	unsigned int n = x.length();
	Rcpp::NumericVector out(n);

	for (unsigned int i = 0; i < n; i++) {
		out(i) = d_target(x(i), kappa, d, log, normalize);
	}

	return out;
}

// [[Rcpp::export]]
inline Rcpp::NumericVector p_target(const Rcpp::NumericVector& q, double kappa,
	double d, bool lower = true, bool log = false)
{
	unsigned int n = q.length();
	Rcpp::NumericVector out(n);

	for (unsigned int i = 0; i < n; i++) {
		out(i) = p_target(q(i), kappa, d, lower, log);
	}

	return out;
}

#endif
