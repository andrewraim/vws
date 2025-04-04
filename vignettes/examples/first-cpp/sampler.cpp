// [[Rcpp::depends(vws, fntl)]]
#include "vws.h"
#include "UniformHelper.h"

// [[Rcpp::export]]
Rcpp::List sample(unsigned int n, double kappa, double d, double lo,
    double hi, unsigned int N)
{
    vws::rejection_args args;
    args.max_rejects = 1000;
    args.report_period = 100;

    const vws::uv_weight_function& w =
    [&](double x, bool log = true) -> double {
        double out = R_NegInf;
        if (lo < x && x <= hi) {
            out = 0.5 * (d - 3) * std::log1p(-std::pow(x, 2)) + kappa*x;
        }
        return log ? out : std::exp(out);
    };

    UniformHelper helper(-1, 1);
    vws::UnivariateConstRegion supp(lo, hi, w, helper);
    vws::FMMProposal<double, vws::UnivariateConstRegion> h({ supp });

    h.adapt(N - 1);
    const auto& out = vws::rejection(h, n, args);

    return Rcpp::List::create(
        Rcpp::Named("draws") = out.draws,
        Rcpp::Named("rejects") = out.rejects
    );
}
