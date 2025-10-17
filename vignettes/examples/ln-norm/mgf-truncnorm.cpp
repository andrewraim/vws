// [[Rcpp::depends(vws, fntl)]]
#include "mgf-truncnorm.h"
#include "vws.h"

double mgf_truncnorm(double s, double a, double b, double z, double lambda,
	bool log)
{
	// If we are truncating way into the upper tail of a distribution, working
	// with the complement of the CDF helps to retain precision. Otherwise,
	// work with the CDF function.

	double lambda2 = lambda * lambda;
	double mean_j = z + s*lambda2;

	double lpa = R::pnorm(a, mean_j, lambda, true, true);
	double lpb = R::pnorm(b, mean_j, lambda, true, true);
	double clpa = R::pnorm(a, mean_j, lambda, false, true);
	double clpb = R::pnorm(b, mean_j, lambda, false, true);
	double lp_num = std::max(
       	vws::log_sub2_exp(clpa, clpb),
       	vws::log_sub2_exp(lpb, lpa)
   	);

	lpa = R::pnorm(a, z, lambda, true, true);
	lpb = R::pnorm(b, z, lambda, true, true);
	clpa = R::pnorm(a, z, lambda, false, true);
	clpb = R::pnorm(b, z, lambda, false, true);
	double lp_den = std::max(
		vws::log_sub2_exp(clpa, clpb),
		vws::log_sub2_exp(lpb, lpa)
	);

	double out = lp_num - lp_den + s*z + s*s * lambda2 / 2;
	return log ? out : std::exp(out);
}
