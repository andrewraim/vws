// [[Rcpp::depends(vws, fntl)]]
#include "vws.h"

// [[Rcpp::export]]
Rcpp::List sample2(unsigned int n, double a, double nu, unsigned int N,
	unsigned int max_rejects = 10000, unsigned int report = 1000)
{
    vws::rejection_args args;
    args.max_rejects = max_rejects;
    args.report = report;

    double lambda = a*a / 4;

    const vws::uv_weight_function& w =
    [&](double x, bool log = true) {
        double out = -lgamma(x + nu + 1);
        return log ? out : std::exp(out);
    };

    fntl::density df = [&](double x, bool log = false) {
        return R::dpois(x, lambda, log);
    };
    fntl::cdf pf = [&](double q, bool lower = true, bool log = false) {
        return R::ppois(q, lambda, lower, log);
    };
    fntl::quantile qf = [&](double p, bool lower = true, bool log = false) {
        return R::qpois(p, lambda, lower, log);
    };

	// Note that w is a decreasing function in x. We can get maximize it if
	// we evaluate w at the smallest integer in the region, and minimize it
	// at the largest integer.

    const vws::optimizer& maxopt = [](const vws::uv_weight_function& w,
    	double lo, double hi, bool log)
    {
    	if (lo < 0 && hi < 0) { Rcpp::stop("Did not code this case"); }
    	double out = (lo < 0) ? w(0, true) : w(std::ceil(lo), true);
    	return log ? out : std::exp(out);
    };

    const vws::optimizer& minopt = [](const vws::uv_weight_function& w,
    	double lo, double hi, bool log)
    {
    	if (lo < 0 && hi < 0) { Rcpp::stop("Did not code this case"); }
    	double out = std::isinf(hi) ? w(hi, true) : w(std::floor(hi), true);
    	return log ? out : std::exp(out);
    };

    vws::UnivariateHelper helper(df, pf, qf);
    vws::IntConstRegion supp(R_NegInf, R_PosInf, w, helper, maxopt, minopt);
    vws::FMMProposal<double, vws::IntConstRegion> h({ supp });

    auto lbdd = h.refine(N - 1);
    auto out = vws::rejection(h, n, args);

    return Rcpp::List::create(
        Rcpp::Named("draws") = out.draws,
        Rcpp::Named("rejects") = out.rejects,
        Rcpp::Named("lbdd") = lbdd
    );
}
