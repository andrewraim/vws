#include "rcpp-interface.h"
#include "vws.h"

Rcpp::NumericVector rect(const Rcpp::NumericVector& z,
	const Rcpp::NumericVector& a, const Rcpp::NumericVector& b)
{
	return vws::rect(z, a, b);
}

Rcpp::NumericVector inv_rect(const Rcpp::NumericVector& x,
	const Rcpp::NumericVector& a, const Rcpp::NumericVector& b)
{
	return vws::inv_rect(x, a, b);
}
