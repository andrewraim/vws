// [[Rcpp::depends(vws, fntl)]]
#include "vws.h"
#include "LinearVWSRegion.h"

// [[Rcpp::export]]
Rcpp::List r_bessel_v3(unsigned int n, double lambda, double nu, unsigned int N,
	double tol = 0, unsigned int max_rejects = 10000, unsigned int report = 10000)
{
	vws::rejection_args args;
	args.max_rejects = max_rejects;
	args.report = report;

	LinearVWSRegion supp(0, R_PosInf, lambda, nu);
	vws::FMMProposal<double, LinearVWSRegion> h(supp);

	auto lbdd = h.refine(N - 2, tol);
	auto out = vws::rejection(h, n, args);

	return Rcpp::List::create(
		Rcpp::Named("draws") = out.draws,
		Rcpp::Named("rejects") = out.rejects,
		Rcpp::Named("lbdd") = lbdd
	);
}
