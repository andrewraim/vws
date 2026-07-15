#ifndef LN_NORM_PROPOSAL_H
#define LN_NORM_PROPOSAL_H

// [[Rcpp::depends(vws, fntl)]]
#include "vws.h"

/*
* Define a subclass of fmm_proposal that we can expose to R via Modules
*/
class ln_norm_proposal : public vws::fmm_proposal<double, vws::real_const_region>
{
public:
	ln_norm_proposal(
		double z,
		double lambda,
		double xbeta,
		double sigma)
	: vws::fmm_proposal<double, vws::real_const_region>(supp(z, lambda, xbeta, sigma))
	{
	}

	void update(double xbeta, double sigma);

	vws::rejection_result<double> draw(unsigned int max_rejects);

	vws::rejection_result<double> draw_tune(double tol_suff, double tol_merge,
		unsigned int max_rejects);

	unsigned int size() const {
		return vws::fmm_proposal<double, vws::real_const_region>::size();
	}

private:
	vws::real_const_region supp(double z, double lambda, double xbeta, double sigma);
};

RCPP_MODULE(vws_module) {
	Rcpp::class_<ln_norm_proposal>("ln_norm_proposal")
	.constructor<double,double,double,double>()
	.method("update", &ln_norm_proposal::update)
	.method("draw", &ln_norm_proposal::draw)
	.method("draw_tune", &ln_norm_proposal::draw_tune)
	.method("size", &ln_norm_proposal::size)
	;
}

/*
* Implementation of member functions is below
*/

inline vws::rejection_result<double>
ln_norm_proposal::draw(unsigned int max_rejects)
{
    vws::rejection_args args;
    args.max_rejects = max_rejects;
    args.report = 1e6;

	return vws::rejection(*this, 1, args);
}

inline vws::rejection_result<double>
ln_norm_proposal::draw_tune(double tol_suff, double tol_merge, unsigned int max_rejects)
{
	vws::rejection_args args;
    args.max_rejects = max_rejects;
    args.report = 1e6;
	args.tol_suff = tol_suff;
	args.tol_merge = tol_merge;

	return vws::rejection_tune(*this, 1, args);
}

inline void ln_norm_proposal::update(double xbeta, double sigma)
{
	const vws::dfdb& w = [=](double y, bool log = true) -> double {
		double out = (y > 0) ? R::dlnorm(y, xbeta, sigma, true) : R_NegInf;
		return log ? out : std::exp(out);
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

	// Update weight function in each region in the proposal
	std::set<vws::real_const_region>::iterator itr = _regions.begin();
	for (; itr != _regions.end(); ++itr) {
		vws::real_const_region& reg = const_cast<vws::real_const_region&>(*itr);
		reg.set_w(w);
		reg.set_maxopt(maxopt);
		reg.set_minopt(minopt);
		reg.init();
	}
}

inline vws::real_const_region ln_norm_proposal::supp(
	double z,
	double lambda,
	double xbeta,
	double sigma)
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
	vws::real_const_region out(0, R_PosInf, w, helper, maxopt, minopt);
	return out;
}

#endif
