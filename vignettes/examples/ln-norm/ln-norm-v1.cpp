// [[Rcpp::depends(vws, fntl)]]
#include "vws.h"

// [[Rcpp::export]]
Rcpp::List r_ln_norm_v1(unsigned int n, double z, double mu,
    double sigma, double lambda, unsigned int N = 10,
    unsigned int max_rejects = 10000, unsigned int report = 1000)
{
    vws::rejection_args args;
    args.max_rejects = max_rejects;
    args.report = report;
    args.action = vws::error_action::STOP;

    const vws::weight_dfd& w =
    [&](double x, bool log = true) {
        double out = R_NegInf;
        if (x > 0) {
            out = -std::log(x) - std::pow(std::log(x) - mu, 2) / (2*sigma*sigma);
        }
        return log ? out : std::exp(out);
    };

    fntl::density df = [&](double x, bool log = false) {
        return R::dnorm(x, z, lambda, log);
    };
    fntl::cdf pf = [&](double q, bool lower = true, bool log = false) {
        return R::pnorm(q, z, lambda, lower, log);
    };
    fntl::quantile qf = [&](double p, bool lower = true, bool log = false) {
        return R::qnorm(p, z, lambda, lower, log);
    };

    vws::UnivariateHelper helper(df, pf, qf);
    vws::RealConstRegion supp(0, R_PosInf, w, helper);
    vws::FMMProposal<double, vws::RealConstRegion> h({ supp });

    h.refine(N - 1);
    h.print(5);

    const vws::rejection_result<double>& out = vws::rejection(h, n, args);

    return Rcpp::List::create(
        Rcpp::Named("draws") = out.draws,
        Rcpp::Named("rejects") = out.rejects
    );
}
