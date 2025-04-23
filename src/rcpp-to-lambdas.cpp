#include "rcpp-to-lambdas.h"

RejectionLambdas rcpp_to_lambdas(const Rcpp::Function& w,
	const Rcpp::Function& d_base, const Rcpp::Function& p_base,
	const Rcpp::Function& q_base)
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

	RejectionLambdas out;
	out.w = w0;
	out.d = d0;
	out.p = p0;
	out.q = q0;
	return out;
}

RejectionLambdas rcpp_to_lambdas(const Rcpp::Function& w,
	const Rcpp::Function& d_base, const Rcpp::Function& p_base,
	const Rcpp::Function& q_base, const Rcpp::Function& maxopt,
	const Rcpp::Function& minopt)
{
	RejectionLambdas out = rcpp_to_lambdas(w, d_base, p_base, q_base);

	const vws::optimizer& maxopt0 = [&](const vws::uv_weight_function& ww,
		double lo, double hi, bool log) -> double
	{
		// Note: here we take w from the captured w rather than the argument
		// `ww`. This helps us avoid having to create an Rcpp::Function that
		// calls `ww`.
		const Rcpp::NumericVector& lo0 = Rcpp::NumericVector::create(lo);
		const Rcpp::NumericVector& hi0 = Rcpp::NumericVector::create(hi);
		const Rcpp::LogicalVector& log0 = Rcpp::LogicalVector::create(log);
		const Rcpp::NumericVector& out = maxopt(w, lo0, hi0, log0);
		return out(0);
	};

	const vws::optimizer& minopt0 = [&](const vws::uv_weight_function& ww,
		double lo, double hi, bool log) -> double
	{
		// Note: here we take w from the captured w rather than the argument
		// `ww`. This helps us avoid having to create an Rcpp::Function that
		// calls `ww`.
		const Rcpp::NumericVector& lo0 = Rcpp::NumericVector::create(lo);
		const Rcpp::NumericVector& hi0 = Rcpp::NumericVector::create(hi);
		const Rcpp::LogicalVector& log0 = Rcpp::LogicalVector::create(log);
		const Rcpp::NumericVector& out = minopt(w, lo0, hi0, log0);
		return out(0);
	};


	out.maxopt = maxopt0;
	out.minopt = minopt0;
	return out;
}
