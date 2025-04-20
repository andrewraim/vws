#include "rcpp-rejection-int.h"
#include "vws.h"
#include "rcpp-to-lambdas.h"

Rcpp::List rejection_int_rcpp(unsigned int n, double lo, double hi,
	const Rcpp::Function& w, const Rcpp::Function& d_base,
	const Rcpp::Function& p_base, const Rcpp::Function& q_base,
	const Rcpp::Function& s_base, unsigned int N, double tol,
	const Rcpp::List& control)
{
	const RejectionLambdas& lam = rcpp_to_lambdas(w, d_base, p_base, q_base, s_base);
	vws::UnivariateHelper helper(lam.d, lam.p, lam.q, lam.s);

	vws::IntConstRegion supp(lo, hi, lam.w, helper);
	vws::FMMProposal<double, vws::IntConstRegion> h({ supp });

	vws::rejection_args args(control);
	const auto& adapt_out = h.adapt(N - 1, tol);
	const auto& out = vws::rejection(h, n, args);

	return Rcpp::List::create(
		Rcpp::Named("draws") = out.draws,
		Rcpp::Named("rejects") = out.rejects,
		Rcpp::Named("lbdd") = adapt_out
	);

	return Rcpp::List::create();
}
