// [[Rcpp::depends(vws, fntl)]]
#include "vws.h"

// [[Rcpp::export]]
Rcpp::List sample2(unsigned int n, double kappa, double d, unsigned int N,
    double tol = 0, unsigned int max_rejects = 10000, unsigned int report = 1000)
{
    vws::rejection_args args;
    args.max_rejects = max_rejects;
    args.report = report;

    const vws::uv_weight_function& w =
    [&](double x, bool log = true) {
        double out = 0.5 * (d - 3) * std::log1p(-std::pow(x, 2)) + kappa*x;
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

    const vws::optimizer& maxopt = [&](const vws::uv_weight_function& w,
        double lo, double hi, bool log)
    {
        double A = 1, B = (d-3) / kappa, C = -1;
        double x_root1 = ( -B + std::sqrt(B*B - 4*A*C) ) / (2*A);
        double x_root2 = ( -B - std::sqrt(B*B - 4*A*C) ) / (2*A);
        double x_root = -1 < x_root1 && x_root1 <= 1 ? x_root1 : x_root2;
        if (-1 >= x_root || x_root > 1) { Rcpp::stop("Invalid root"); }

        double out;
        if (x_root < lo) {
            out = w(lo, true);
        } else if (x_root > hi) {
            out = w(hi, true);
        } else {
            out = w(x_root, true);
        }

        return out;
    };

    const vws::optimizer& minopt = [&](const vws::uv_weight_function& w,
        double lo, double hi, bool log)
    {
        return w(std::min(lo, hi), true);
    };

    vws::UnivariateHelper helper(df, pf, qf);
    vws::RealConstRegion supp(-1, 1, w, helper, maxopt, minopt);
    vws::FMMProposal<double, vws::RealConstRegion> h({ supp });

    auto lbdd = h.refine(N - 1);
    auto out = vws::rejection(h, n, args);

    return Rcpp::List::create(
        Rcpp::Named("draws") = out.draws,
        Rcpp::Named("rejects") = out.rejects,
        Rcpp::Named("lbdd") = lbdd
    );
}
