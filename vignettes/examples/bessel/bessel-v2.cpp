// [[Rcpp::depends(vws, fntl)]]
#include "vws.h"

// [[Rcpp::export]]
Rcpp::List r_bessel_v2(unsigned int n, double lambda, double nu, unsigned int N,
	double tol = 0, unsigned int max_rejects = 10000, unsigned int report = 10000,
	const Rcpp::NumericVector& x = Rcpp::NumericVector::create())
{
	vws::rejection_args args;
	args.max_rejects = max_rejects;
	args.report = report;

	double mean = lambda * lambda / 4;

	const vws::dfdb& w =
	[&](double x, bool log = true) {
		double out = -std::lgamma(x + nu + 1);
		return log ? out : std::exp(out);
	};

	fntl::density df = [&](double x, bool log = false) {
		return R::dpois(x, mean, log);
	};
	fntl::cdf pf = [&](double q, bool lower = true, bool log = false) {
		return R::ppois(q, mean, lower, log);
	};
	fntl::quantile qf = [&](double p, bool lower = true, bool log = false) {
		return R::qpois(p, mean, lower, log);
	};

	// Note that w is a decreasing function in x. We can get maximize it if
	// we evaluate w at the smallest integer in the region, and minimize it
	// at the largest integer. We can also ensure that only integer values
	// are considered, which is more awkward with numerical optimization.

	const vws::optimizer& maxopt = [](const vws::dfdb& w,
		double lo, double hi, bool log)
	{
		if (lo < 0 && hi < 0) { Rcpp::stop("Did not code this case"); }
		double x = (lo < 0) ? 0 : std::floor(lo) + 1;
		double out = w(x, true);
		return log ? out : std::exp(out);
	};

	const vws::optimizer& minopt = [](const vws::dfdb& w,
		double lo, double hi, bool log)
	{
		if (lo < 0 && hi < 0) { Rcpp::stop("Did not code this case"); }
		double x = std::isinf(hi) ? hi : std::floor(hi);
		double out = w(x, true);
		return log ? out : std::exp(out);
	};

	vws::UnivariateHelper helper(df, pf, qf);
	vws::IntConstRegion supp(-0.1, R_PosInf, w, helper, maxopt, minopt);
	vws::FMMProposal<double, vws::IntConstRegion> h(supp);

	auto lbdd = h.refine(N - 1, tol);
	auto out = vws::rejection(h, n, args);

	/*
	* Evaluate the majorized weight function, minorized weight function, and
	* proposal on the given values of x. The caller can plot these.
	*/
	Rcpp::NumericVector hx(x.size());
	for (unsigned int i = 0; i < x.size(); i++) {
		hx(i) = h.d(x(i));
	}

	Rcpp::NumericVector lower(N);
	Rcpp::NumericVector upper(N);
	Rcpp::NumericVector lwmin(N);
	Rcpp::NumericVector lwmax(N);

	std::set<vws::IntConstRegion>::const_iterator itr = h.regions_begin();
	unsigned int i = 0;
	for (; itr != h.regions_end(); itr++) {
		lower(i) = itr->lower();
		upper(i) = itr->upper();
		lwmin(i) = itr->w_minor(itr->midpoint());
		lwmax(i) = itr->w_major(itr->midpoint());
		i++;
	}

	const Rcpp::DataFrame& df_weight = Rcpp::DataFrame::create(
		Rcpp::Named("lo") = lower,
		Rcpp::Named("hi") = upper,
		Rcpp::Named("lwmin") = lwmin,
		Rcpp::Named("lwmax") = lwmax
	);

	return Rcpp::List::create(
		Rcpp::Named("draws") = out.draws,
		Rcpp::Named("rejects") = out.rejects,
		Rcpp::Named("lbdd") = lbdd,
		Rcpp::Named("hx") = hx,
		Rcpp::Named("df_weight") = df_weight
	);
}
