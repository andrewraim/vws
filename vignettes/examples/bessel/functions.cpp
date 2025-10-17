// [[Rcpp::depends(vws, fntl)]]
#include "functions.h"
#include "vws.h"

double mgf_truncpois(double s, double lambda, double a, double b, bool log)
{
	// If we are truncating way into the upper tail of a distribution, working
	// with the complement of the CDF helps to retain precision. Otherwise,
	// work with the CDF function.

	double lo = std::max(a, 0.0);
	double hi = std::max(b, 0.0);
	double lambda2 = lambda * lambda;
	double mean = std::exp(s) * lambda2 / 4;

	double lpa = R::ppois(lo, mean, true, true);
	double lpb = R::ppois(hi, mean, true, true);
	double clpa = R::ppois(lo, mean, false, true);
	double clpb = R::ppois(hi, mean, false, true);
	double lp = vws::log_sub2_exp(lpb, lpa);
	double clp = vws::log_sub2_exp(clpa, clpb);
	double lm = std::max(lp, clp);

	// double lg1 = incgamma(std::floor(b) + 1, mean, false, true);
	// double lg2 = incgamma(std::floor(b) + 1, 0, false, true);
	// double lg3 = incgamma(std::ceil(a), mean, false, true);
	// double lg4 = incgamma(std::ceil(a), 0, false, true);
	// double lm = vws::log_sub2_exp(lg1 - lg2, lg3 - lg4);

	// double out = lambda2 / 4 * std::expm1(s) + lm;
	double out = -lambda2 / 4 + mean + lm;
	return log ? out : std::exp(out);
}

Rcpp::NumericVector incgamma(double a, const Rcpp::NumericVector& x, bool lower,
	bool log)
{
	unsigned int n = x.size();
	Rcpp::NumericVector out(n);

	for (unsigned int i = 0; i < n; i++) {
		out(i) = incgamma(a, x(i), lower, log);
	}

	return out;
}

double incgamma(double a, double x, bool lower, bool log)
{
	double out = std::lgamma(a) + R::pgamma(x, a, 1, lower, true);
	return log ? out : std::exp(out);
}
