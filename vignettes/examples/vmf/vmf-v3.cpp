// [[Rcpp::depends(vws, fntl)]]
#include "vws.h"
#include "LinearVWSRegion.h"

// [[Rcpp::export]]
Rcpp::List r_vmf_pre_v3(unsigned int n, double kappa, double d,
	unsigned int N, double tol, unsigned int max_rejects, unsigned int report)
{
	vws::rejection_args args;
	args.max_rejects = max_rejects;
	args.report = report;
	args.action = vws::error_action::STOP;

	LinearVWSRegion supp(-1, 1, kappa, d);
	vws::FMMProposal<double, LinearVWSRegion> h({ supp });

	auto lbdd = h.refine(N - 1, tol);
	auto out = vws::rejection(h, n, args);

	return Rcpp::List::create(
		Rcpp::Named("draws") = out.draws,
		Rcpp::Named("rejects") = out.rejects,
		Rcpp::Named("lbdd") = lbdd
	);
}
