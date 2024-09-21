#ifndef VWS_GUMBEL_H
#define VWS_GUMBEL_H

#include <Rcpp.h>

namespace vws {

/*
* Non-vectorized distribution functions.
*/

inline double q_gumbel(double p, double mu = 0, double sigma = 1,
	bool lower = true, bool log = false)
{
	double lp0 = log ? p : std::log(p);
	double lp = lower ? lp0 : std::log1p(-std::exp(lp0));
	return mu - sigma * std::log(-lp);
}

inline double r_gumbel(double mu = 0, double sigma = 1)
{
	double u = R::runif(0, 1);
	return q_gumbel(u, mu, sigma);
}

inline double d_gumbel(double x, double mu = 0,
	double sigma = 1, bool log = false)
{
	double z = (x - mu) / sigma;
	double out = -std::log(sigma) - (z + std::exp(-z));
	return log ? out : std::exp(out);
}

inline double p_gumbel(double q, double mu = 0,
	double sigma = 1, bool lower = true, bool log = false)
{
	double z = (q - mu) / sigma;
	double out0 = -std::exp(-z);
	double out = lower ? out0 : std::log1p(-std::exp(out0));
	return log ? out : std::exp(out);
}

/*
* Vectorized distribution functions.
*/

inline Rcpp::NumericVector r_gumbel(unsigned int n, double mu = 0, double sigma = 1)
{
	const Rcpp::NumericVector& u = Rcpp::runif(n);
	Rcpp::NumericVector out(n);

	for (unsigned int i = 0; i < n; i++) {
		out(i) = q_gumbel(u(i), mu, sigma);
	}

	return out;
}

inline Rcpp::NumericVector d_gumbel(const Rcpp::NumericVector& x, double mu = 0,
	double sigma = 1, bool log = false)
{
	unsigned int n = x.size();
	Rcpp::NumericVector out(n);

	for (unsigned int i = 0; i < n; i++) {
		out(i) = d_gumbel(x(i), mu, sigma, log);
	}

	return out;
}

inline Rcpp::NumericVector p_gumbel(const Rcpp::NumericVector& q, double mu = 0,
	double sigma = 1, bool lower = true, bool log = false)
{
	unsigned int n = q.size();
	Rcpp::NumericVector out(n);

	for (unsigned int i = 0; i < n; i++) {
		out(i) = p_gumbel(q(i), mu, sigma, lower, log);
	}

	return out;
}

inline Rcpp::NumericVector q_gumbel(const Rcpp::NumericVector& p, double mu = 0,
	double sigma = 1, bool lower = true, bool log = false)
{
	unsigned int n = p.size();
	Rcpp::NumericVector out(n);

	for (unsigned int i = 0; i < n; i++) {
		out(i) = q_gumbel(p(i), mu, sigma, lower, log);
	}

	return out;
}

}

#endif
