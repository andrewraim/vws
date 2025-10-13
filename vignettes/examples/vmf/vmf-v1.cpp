// [[Rcpp::depends(vws, fntl)]]                             // <1>
#include "vws.h"                                            // <2>
#include "texp.h"                                           // <3>

// [[Rcpp::export]]                                         // <4>
Rcpp::List r_vmf_pre_v1(unsigned int n, double kappa, double d,
	unsigned int N, double tol = 0, unsigned int max_rejects = 10000,
	unsigned int report = 10000)
{
    vws::rejection_args args;                               // <5>
    args.max_rejects = max_rejects;
    args.report = report;

    const vws::weight_dfd& w =                              // <6>
    [&](double x, bool log = true) {
        double out = R_NegInf;
        if (std::fabs(x) < 1){
            out = 0.5 * (d - 3) * std::log1p(-std::pow(x, 2));
        }
        return log ? out : std::exp(out);
    };

    fntl::density df = [&](double x, bool log = false) {    // <7>
        return d_texp(x, kappa, -1, 1, log);
    };
    fntl::cdf pf = [&](double q, bool lower = true, bool log = false) {
        return p_texp(q, kappa, -1, 1, lower, log);
    };
    fntl::quantile qf = [&](double p, bool lower = true, bool log = false) {
        return q_texp(p, kappa, -1, 1, lower, log);
    };

    vws::UnivariateHelper helper(df, pf, qf);               // <8>
    vws::RealConstRegion supp(-1, 1, w, helper);            // <9>
    vws::FMMProposal<double, vws::RealConstRegion> h(supp); // <10>

    auto lbdd = h.refine(N - 1, tol);                       // <11>
    auto out = vws::rejection(h, n, args);                  // <12>

    return Rcpp::List::create(                              // <13>
        Rcpp::Named("draws") = out.draws,
        Rcpp::Named("rejects") = out.rejects,
        Rcpp::Named("lbdd") = lbdd
    );
}
