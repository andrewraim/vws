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
		double y,
		double lambda,
		double xbeta,
		double sigma)
	: vws::fmm_proposal<double, vws::real_const_region>(supp(y, lambda, xbeta, sigma))
	{
	}

	void update(double xbeta, double sigma);

private:
	vws::real_const_region supp(double y, double lambda, double xbeta, double sigma);
};

/*
* Implementation of member functions is below
*/

inline void ln_norm_proposal::update(double xbeta, double sigma)
{
	const vws::dfdb& w = [=](double mu, bool log = true) -> double {
		double out = (mu > 0) ? R::dlnorm(mu, xbeta, sigma, true) : R_NegInf;
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
	double y,
	double lambda,
	double xbeta,
	double sigma)
{
	const vws::dfdb& w = [=](double mu, bool log = true) -> double {
		double out = (mu > 0) ? R::dlnorm(mu, xbeta, sigma, true) : R_NegInf;
		return log ? out : std::exp(out);
	};

	const fntl::density& df = [=](double mu, bool log = false) {
		return R::dnorm(mu, y, lambda, log);
	};

	const fntl::cdf& pf = [=](double q, bool lower = true, bool log = false) {
		return R::pnorm(q, y, lambda, lower, log);
	};

	const fntl::quantile& qf = [=](double p, bool lower = true, bool log = false) {
		return R::qnorm(p, y, lambda, lower, log);
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
