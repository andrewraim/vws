// [[Rcpp::depends(vws, fntl)]]
#include "vws.h"

// [[Rcpp::export]]
Rcpp::List sample(unsigned int n, double kappa, double d, unsigned int N)
{
    vws::rejection_args args;
    args.max_rejects = 1000;
    args.report = 100;

    const vws::uv_weight_function& w =
    [&](double x, bool log = true) -> double {
        double out = R_NegInf;
        if (-1 < x && x <= 1) {
            out = 0.5 * (d - 3) * std::log1p(-std::pow(x, 2)) + kappa*x;
        }
        return log ? out : std::exp(out);
    };

    fntl::density df = [](double x, bool log = false) {
        return R::dunif(x, -1, 1, log);
    };
    fntl::cdf pf = [](double q, bool lower = true, bool log = false) {
        return R::punif(q, -1, 1, lower, log);
    };
    fntl::quantile qf = [](double p, bool lower = true, bool log = false) {
        return R::qunif(p, -1, 1, lower, log);
    };
    vws::supp sf = [](double x) {
        return -1 <= x && x <= 1;
    };

    vws::UnivariateHelper helper(df, pf, qf, sf);
    vws::UnivariateConstRegion supp(-1, 1, w, helper);
    vws::FMMProposal<double, vws::UnivariateConstRegion> h({ supp });

    h.adapt(N - 1);
    const auto& out = vws::rejection(h, n, args);

    return Rcpp::List::create(
        Rcpp::Named("draws") = out.draws,
        Rcpp::Named("rejects") = out.rejects
    );
}
