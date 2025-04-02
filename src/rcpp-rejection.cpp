#include "rcpp-rejection.h"
#include "vws.h"

Rcpp::List rejection_rcpp(unsigned int n, double lo, double hi,
	const Rcpp::Function& w, const Rcpp::Function& d_base,
	const Rcpp::Function& p_base, const Rcpp::Function& q_base,
	const Rcpp::Function& s_base, unsigned int N, double tol,
	const Rcpp::List& control)
{
	const vws::uv_weight_function& w0 =
	[&](double x, bool log = true) -> double {
		const Rcpp::NumericVector& x0 = Rcpp::NumericVector::create(x);
		const Rcpp::LogicalVector& log0 = Rcpp::LogicalVector::create(log);
		const Rcpp::NumericVector& out = w(x0, log0);
		return out(0);
	};

	const fntl::density& d0 =
	[&](double x, bool log = true) -> double {
		const Rcpp::NumericVector& x0 = Rcpp::NumericVector::create(x);
		const Rcpp::LogicalVector& log0 = Rcpp::LogicalVector::create(log);
		const Rcpp::NumericVector& out = d_base(x0, log0);
		return out(0);
	};

	const fntl::cdf& p0 =
	[&](double x, bool lower = true, bool log = true) -> double {
		const Rcpp::NumericVector& x0 = Rcpp::NumericVector::create(x);
		const Rcpp::LogicalVector& lower0 = Rcpp::LogicalVector::create(lower);
		const Rcpp::LogicalVector& log0 = Rcpp::LogicalVector::create(log);
		const Rcpp::NumericVector& out = p_base(x0, lower0, log0);
		return out(0);
	};

	const fntl::quantile& q0 =
	[&](double x, bool lower = true, bool log = true) -> double {
		const Rcpp::NumericVector& x0 = Rcpp::NumericVector::create(x);
		const Rcpp::LogicalVector& lower0 = Rcpp::LogicalVector::create(lower);
		const Rcpp::LogicalVector& log0 = Rcpp::LogicalVector::create(log);
		const Rcpp::NumericVector& out = q_base(x0, lower0, log0);
		return out(0);
	};

	const vws::supp& s0 =
	[&](double x) -> double {
		const Rcpp::NumericVector& x0 = Rcpp::NumericVector::create(x);
		const Rcpp::NumericVector& out = s_base(x0);
		return out(0);
	};

	vws::LambdaHelper helper(lo, hi, d0, p0, q0, s0);
	vws::UnivariateConstRegion supp(lo, hi, w0, helper);
	vws::FMMProposal<double, vws::UnivariateConstRegion> h({ supp });

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
