// [[Rcpp::depends(vws, fntl)]]
#include "mgf-truncpois.h"
#include "vws.h"

double mgf_truncpois(double s, double a, double b, double lambda,
	bool log)
{
	// If we are truncating way into the upper tail of a distribution, working
	// with the complement of the CDF helps to retain precision. Otherwise,
	// work with the CDF function.

	double lambda2 = lambda * lambda;
	double mean = std::exp(s) * lambda2 / 4;

	/*
	double lpa = R::ppois(a, mean, true, true);
	double lpb = R::ppois(b, mean, true, true);
	double clpa = R::ppois(a, mean, false, true);
	double clpb = R::ppois(b, mean, false, true);
	double lp_num = std::max(
       	vws::log_sub2_exp(clpa, clpb),
       	vws::log_sub2_exp(lpb, lpa)
   	);

	lpa = R::ppois(a, mean, true, true);
	lpb = R::ppois(b, mean, true, true);
	clpa = R::ppois(a, mean, false, true);
	clpb = R::ppois(b, mean, false, true);
	double lp_den = std::max(
		vws::log_sub2_exp(clpa, clpb),
		vws::log_sub2_exp(lpb, lpa)
	);
	*/

	double lg1 = incgamma(std::floor(b) + 1, mean, false, true);
	double lg2 = incgamma(std::floor(b) + 1, 0, false, true);
	double lg3 = incgamma(std::ceil(a), mean, false, true);
	double lg4 = incgamma(std::ceil(a), 0, false, true);
	double lm = vws::log_sub2_exp(lg1 - lg2, lg3 - lg4);

	double out = lambda2 / 4 * std::expm1(s) + lm;
	return log ? out : std::exp(out);
}

double incgamma(double a, double x, bool lower, bool log)
{
	double out = std::lgamma(a) + R::pgamma(x, a, 1, !lower, true);
	return log ? out : std::exp(out);
}
