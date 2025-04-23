#include "rcpp-rejection-int.h"
#include "vws.h"
#include "rcpp-to-lambdas.h"

Rcpp::List rejection_int_rcpp(unsigned int n, double lo, double hi,
	const Rcpp::Function& w, const Rcpp::Function& d_base,
	const Rcpp::Function& p_base, const Rcpp::Function& q_base,
	unsigned int N, double tol, const Rcpp::List& control)
{
	const RejectionLambdas& lam = rcpp_to_lambdas(w, d_base, p_base, q_base);
	vws::UnivariateHelper helper(lam.d, lam.p, lam.q);
	vws::IntConstRegion supp(lo, hi, lam.w, helper);
	vws::FMMProposal<double, vws::IntConstRegion> h({ supp });

	vws::rejection_args args(control);
	const auto& adapt_out = h.adapt(N - 1, tol);
	const auto& rejection_out = vws::rejection(h, n, args);

	return Rcpp::List::create(
		Rcpp::Named("draws") = rejection_out.draws,
		Rcpp::Named("rejects") = rejection_out.rejects,
		Rcpp::Named("lbdd") = adapt_out
	);
}

Rcpp::List rejection_int_opt_rcpp(unsigned int n, double lo, double hi,
	const Rcpp::Function& w, const Rcpp::Function& d_base,
	const Rcpp::Function& p_base, const Rcpp::Function& q_base,
	const Rcpp::Function& maxopt, const Rcpp::Function& minopt,
	unsigned int N, double tol, const Rcpp::List& control)
{
	const RejectionLambdas& lam = rcpp_to_lambdas(w, d_base, p_base, q_base,
		maxopt, minopt);
	vws::UnivariateHelper helper(lam.d, lam.p, lam.q);
	vws::IntConstRegion supp(lo, hi, lam.w, helper, lam.maxopt, lam.minopt);
	vws::FMMProposal<double, vws::IntConstRegion> h({ supp });

	vws::rejection_args args(control);
	const auto& adapt_out = h.adapt(N - 1, tol);
	const auto& rejection_out = vws::rejection(h, n, args);

	return Rcpp::List::create(
		Rcpp::Named("draws") = rejection_out.draws,
		Rcpp::Named("rejects") = rejection_out.rejects,
		Rcpp::Named("lbdd") = adapt_out
	);
}
