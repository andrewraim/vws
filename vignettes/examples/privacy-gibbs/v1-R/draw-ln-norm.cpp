// [[Rcpp::depends(vws, fntl)]]
#include "vws.h"

// [[Rcpp::export]]
Rcpp::List draw_ln_norm(
	double z,
	double lambda,
	double xbeta,
	double sigma,
	double tol,
	unsigned int N,
	unsigned int max_rejects)
{
	const vws::dfdb& w = [=](double y, bool log = true) -> double {
		double out = (y > 0) ? R::dlnorm(y, xbeta, sigma, true) : R_NegInf;
		return log ? out : std::exp(out);
	};

	const fntl::density& df = [=](double y, bool log = false) {
		return R::dnorm(y, z, lambda, log);
	};

	const fntl::cdf& pf = [=](double q, bool lower = true, bool log = false) {
		return R::pnorm(q, z, lambda, lower, log);
	};

	const fntl::quantile& qf = [=](double p, bool lower = true, bool log = false) {
		return R::qnorm(p, z, lambda, lower, log);
	};

	// We have a simple closed-form max and min that we can use here.
	double sigma2 = sigma * sigma;
	double mode = exp(xbeta - sigma2);

	const vws::optimizer& maxopt =
	[=](const vws::dfdb& w, double lo, double hi, bool log) {
		double y = (mode > hi) ? hi :
			(mode < lo) ? lo :
			mode;
		double out = w(y, true);
		return log ? out : exp(out);
	};

	const vws::optimizer& minopt =
	[=](const vws::dfdb& w, double lo, double hi, bool log) {
		double lwa = w(lo, true);
		double lwb = w(hi, true);
		double out = std::min(lwa, lwb);
		return log ? out : exp(out);
	};

	vws::univariate_helper helper(df, pf, qf);
	vws::real_const_region supp(0, R_PosInf, w, helper, maxopt, minopt);
	vws::fmm_proposal<double, vws::real_const_region> h({supp});
	h.refine(N - 1, tol);

	vws::rejection_args args;
	args.max_rejects = max_rejects;
	args.report = 1e6;

	const auto& out = vws::rejection(h, 1, args);
	return Rcpp::wrap(out);
}
