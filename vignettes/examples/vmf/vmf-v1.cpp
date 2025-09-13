// [[Rcpp::depends(vws, fntl)]]
#include "vws.h"
#include "texp.h"

// [[Rcpp::export]]
Rcpp::List r_vmf_pre_v1(unsigned int n, double kappa, double d, unsigned int N,
    double tol = 0, unsigned int max_rejects = 10000, unsigned int report = 10000)
{
    vws::rejection_args args;
    args.max_rejects = max_rejects;
    args.report = report;

    const vws::weight_dfd& w =
    [&](double x, bool log = true) {
        double out = R_NegInf;
        if (std::fabs(x) < 1){
            out = 0.5 * (d - 3) * std::log1p(-std::pow(x, 2));
        }
        return log ? out : std::exp(out);
    };

    fntl::density df = [&](double x, bool log = false) {
        return d_texp(x, kappa, -1, 1, log);
    };
    fntl::cdf pf = [&](double q, bool lower = true, bool log = false) {
        return p_texp(q, kappa, -1, 1, lower, log);
    };
    fntl::quantile qf = [&](double p, bool lower = true, bool log = false) {
        return q_texp(p, kappa, -1, 1, lower, log);
    };

    vws::UnivariateHelper helper(df, pf, qf);
    vws::RealConstRegion supp(-1, 1, w, helper);
    vws::FMMProposal<double, vws::RealConstRegion> h({ supp });

    auto lbdd = h.refine(N - 1);
    auto out = vws::rejection(h, n, args);

    return Rcpp::List::create(
        Rcpp::Named("draws") = out.draws,
        Rcpp::Named("rejects") = out.rejects,
        Rcpp::Named("lbdd") = lbdd
    );
}
