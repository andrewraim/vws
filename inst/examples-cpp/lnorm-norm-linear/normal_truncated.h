// [[Rcpp::depends(vws)]]
#ifndef NORMAL_TRUNCATED_H
#define NORMAL_TRUNCATED_H

#include <Rcpp.h>
#include "vws.h"

// [[Rcpp::export]]
double d_norm_trunc(double x, double mu, double sigma, double a, double b, bool log)
{
	double log_num = a < x && x <= b ? R::dnorm(x, mu, sigma, true) : R_NegInf;

	// Compute two ways: using upper tail and lower tail CDF. One of these may
	// be -Inf if x is far into the tail of N(mu, sigma2); in this case, the
	// other one may be finite.

	double log_den1 = vws::log_sub2_exp(
		R::pnorm(b, mu, sigma, true, true),
		R::pnorm(a, mu, sigma, true, true)
	);
	double log_den2 = vws::log_sub2_exp(
		R::pnorm(a, mu, sigma, false, true),
		R::pnorm(b, mu, sigma, false, true)
	);
	double log_den = std::max(log_den1, log_den2);

	double out = log_num - log_den;
	return log ? out : exp(out);
}

// [[Rcpp::export]]
double p_norm_trunc(double x, double mu, double sigma,
	double a, double b, bool lower, bool log)
{
    double out;

	if (a < x && x <= b) {
		double lp_a = R::pnorm(a, mu, sigma, true, true);
		double lp_b = R::pnorm(b, mu, sigma, true, true);
		double lp_x = R::pnorm(x, mu, sigma, true, true);
		double lp_den = vws::log_sub2_exp(lp_b, lp_a);

		double clp_a = R::pnorm(a, mu, sigma, false, true);
		double clp_b = R::pnorm(b, mu, sigma, false, true);
		double clp_x = R::pnorm(x, mu, sigma, false, true);
		double clp_den = vws::log_sub2_exp(clp_a, clp_b);

		out = std::isfinite(lp_den) ?
			vws::log_sub2_exp(lp_x, lp_a) - lp_den :
			vws::log_sub2_exp(clp_a, clp_x) - clp_den;
	} else {
		out = R_NegInf;
	}

	return log ? out : exp(out);
}

// [[Rcpp::export]]
double q_norm_trunc(double p, double mu, double sigma,
	double a, double b, bool lower, bool log)
{
	double lp = log ? p : std::log(p);

	double lp_a = R::pnorm(a, mu, sigma, true, true);
	double lp_b = R::pnorm(b, mu, sigma, true, true);
	double lp_diff = vws::log_sub2_exp(lp_b, lp_a);

	double clp_a = R::pnorm(a, mu, sigma, false, true);
	double clp_b = R::pnorm(b, mu, sigma, false, true);
	double clp_diff = vws::log_sub2_exp(clp_a, clp_b);

	double out;

	if (std::isfinite(lp_diff)) {
		double lp_adj = vws::log_add2_exp(lp_a, lp + vws::log_sub2_exp(lp_b, lp_a));
		out = R::qnorm(lp_adj, mu, sigma, true, true);
	} else {
		double clp_adj = vws::log_sub2_exp(clp_a, lp + vws::log_sub2_exp(clp_a, clp_b));
		out = R::qnorm(clp_adj, mu, sigma, false, true);
	}

    return out;
}

// [[Rcpp::export]]
Rcpp::NumericVector r_norm_trunc(unsigned int n, double mu, double sigma,
	double a, double b)
{
	Rcpp::NumericVector out(n);

	const Rcpp::NumericVector& u = Rcpp::runif(n);
	for (unsigned int i = 0; i < n; i++) {
		out(i) = q_norm_trunc(u(i), mu, sigma, a, b, true, false);
	}

	return out;
}

#endif
