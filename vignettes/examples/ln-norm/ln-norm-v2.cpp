// [[Rcpp::depends(vws, fntl)]]
#include "vws.h"

// [[Rcpp::export]]
Rcpp::List r_ln_norm_v2(unsigned int n, double z, double mu,
    double sigma, double lambda, unsigned int N = 10,
    unsigned int max_rejects = 10000, unsigned int report = 1000)
{
    vws::rejection_args args;
    args.max_rejects = max_rejects;
    args.report = report;
    args.action = vws::error_action::STOP;
    double sigma2 = sigma*sigma;

    const vws::weight_dfd& w =
    [&](double x, bool log = true) {
        double out = R_NegInf;
        if (x > 0) {
            out = -std::log(x) - std::pow(std::log(x) - mu, 2) / (2*sigma2);
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

    const vws::optimizer& maxopt = [&](const vws::weight_dfd& w, double lo,
        double hi, bool log)
    {
        double y_star = exp(mu - sigma2);
        double out;

        if (y_star > hi) {
            out = w(hi, true);
        } else if (y_star < lo) {
            out = w(lo, true);
        } else {
            out = w(y_star, true);
        }

        return log ? out : exp(out);
    };

    const vws::optimizer& minopt = [&](const vws::weight_dfd& w, double lo,
        double hi, bool log)
    {
        double y_star = exp(mu - sigma2);
        double lwa = w(lo, true);
        double lwb = w(hi, true);
        double out = std::min(lwa, lwb);
        return log ? out : exp(out);
    };

    vws::UnivariateHelper helper(df, pf, qf);
    vws::RealConstRegion supp(0, R_PosInf, w, helper, maxopt, minopt);
    vws::FMMProposal<double, vws::RealConstRegion> h({ supp });

    auto lbdd = h.refine(N - 1);
    const vws::rejection_result<double>& out = vws::rejection(h, n, args);

    return Rcpp::List::create(
        Rcpp::Named("draws") = out.draws,
        Rcpp::Named("rejects") = out.rejects,
        Rcpp::Named("lbdd") = lbdd
    );
}
