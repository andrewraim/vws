// [[Rcpp::depends(vws, fntl)]]
#include "vws.h"

// [[Rcpp::export]]
Rcpp::List r_bessel_v1(unsigned int n, double lambda, double nu, unsigned int N,
	double tol = 0, unsigned int max_rejects = 10000, unsigned int report = 10000,
	const Rcpp::NumericVector& x = Rcpp::NumericVector::create())
{
    vws::rejection_args args;
    args.max_rejects = max_rejects;
    args.report = report;

    double mean = lambda * lambda / 4;

    const vws::weight_dfd& w =
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

    vws::UnivariateHelper helper(df, pf, qf);
    vws::IntConstRegion supp(-0.1, R_PosInf, w, helper);
    vws::FMMProposal<double, vws::IntConstRegion> h(supp);

    auto lbdd = h.refine(N - 1, tol);
    auto out = vws::rejection(h, n, args);

    // Evaluate the majorized weight function, minorized weight function, and
    // proposal on the given values of x. The caller can plot these.
	unsigned int m = x.size();
    Rcpp::NumericVector hx(m);
    Rcpp::NumericVector wmajx(m);
    for (unsigned int i = 0; i < m; i++) {
        hx(i) = h.d(x(i));
        wmajx(i) = h.w_major(x(i));
    }

    Rcpp::NumericVector lo(N);
    Rcpp::NumericVector hi(N);
    Rcpp::NumericVector wmin(N);
    Rcpp::NumericVector wmax(N);

    std::set<vws::IntConstRegion>::const_iterator itr = h.regions_begin();
    unsigned int i = 0;
    for (; itr != h.regions_end(); itr++) {
		lo(i) = itr->lower();
    	hi(i) = itr->upper();
    	wmin(i) = itr->w_minor(itr->midpoint());
    	wmax(i) = itr->w_major(itr->midpoint());
    	i++;
    }

    const Rcpp::DataFrame& df_weight = Rcpp::DataFrame::create(
    	Rcpp::Named("lo") = lo,
    	Rcpp::Named("hi") = hi,
    	Rcpp::Named("wmin") = wmin,
    	Rcpp::Named("wmax") = wmax
    );

    return Rcpp::List::create(
        Rcpp::Named("draws") = out.draws,
        Rcpp::Named("rejects") = out.rejects,
        Rcpp::Named("lbdd") = lbdd,
        Rcpp::Named("hx") = hx,
        Rcpp::Named("wmajx") = wmajx,
        Rcpp::Named("df_weight") = df_weight
    );
}
