#include "rcpp-interface.h"
#include "vws.h"
#include "fntl.h"

Rcpp::NumericVector rect_rcpp(const Rcpp::NumericVector& z,
	const Rcpp::NumericVector& a, const Rcpp::NumericVector& b)
{
	return vws::rect(z, a, b);
}

Rcpp::NumericVector inv_rect_rcpp(const Rcpp::NumericVector& x,
	const Rcpp::NumericVector& a, const Rcpp::NumericVector& b)
{
	return vws::inv_rect(x, a, b);
}

Rcpp::List optimize_hybrid_rcpp(const Rcpp::Function& f, double init, double lower,
	double upper, bool maximize, unsigned maxiter)
{
	const fntl::dfd& ff = [&](double x) -> double {
		const Rcpp::NumericVector& out = f(x);
		return out(0);
	};

	const auto& out = vws::optimize_hybrid(ff, init, lower, upper, maximize, maxiter);

	return Rcpp::List::create(
		Rcpp::Named("par") = out.par,
		Rcpp::Named("value") = out.value,
		Rcpp::Named("method") = out.method
	);
}

