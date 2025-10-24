// [[Rcpp::depends(vws, fntl)]]
#include "vws.h"
#include "LinearVWSRegion.h"

// [[Rcpp::export]]
Rcpp::List r_ln_norm_v3(unsigned int n, double z, double mu, double sigma,
    double lambda, double lo, double hi, unsigned int N, double tol = 0,
    unsigned int max_rejects = 10000, unsigned int report = 10000)
{
    vws::rejection_args args;
    args.max_rejects = max_rejects;
    args.report = report;

    // Initially partition at y0 where convexity of the weight function changes
    double y0 = exp(mu - std::pow(sigma, 2) + 1);
    LinearVWSRegion r1(lo, y0, z, mu, sigma, lambda);
    LinearVWSRegion r2(y0, hi, z, mu, sigma, lambda);
    vws::FMMProposal<double, LinearVWSRegion> h({r1, r2});

    auto lbdd = h.refine(N - 2, tol);
    auto out = vws::rejection(h, n, args);

    return Rcpp::List::create(
        Rcpp::Named("draws") = out.draws,
        Rcpp::Named("rejects") = out.rejects,
        Rcpp::Named("lbdd") = lbdd
    );
}
