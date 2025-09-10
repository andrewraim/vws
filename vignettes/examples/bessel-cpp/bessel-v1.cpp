// [[Rcpp::depends(vws, fntl)]]
#include "vws.h"

// [[Rcpp::export]]
Rcpp::List r_bessel_v1(unsigned int n, double a, double nu, unsigned int N,
	unsigned int max_rejects = 10000, unsigned int report = 1000)
{
    vws::rejection_args args;
    args.max_rejects = max_rejects;
    args.report = report;

    double lambda = a*a / 4;

    const vws::weight_dfd& w =
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

    vws::UnivariateHelper helper(df, pf, qf);
    vws::IntConstRegion supp(R_NegInf, R_PosInf, w, helper);
    vws::FMMProposal<double, vws::IntConstRegion> h({ supp });

    auto lbdd = h.refine(N - 1);
    auto out = vws::rejection(h, n, args);

    return Rcpp::List::create(
        Rcpp::Named("draws") = out.draws,
        Rcpp::Named("rejects") = out.rejects,
        Rcpp::Named("lbdd") = lbdd
    );
}
