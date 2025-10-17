#include <Rcpp.h>
#include "vws.h"

Rcpp::NumericVector rcpp_seq(double lo, double hi, unsigned int N, bool endpoints)
{
	std::vector<double> out = vws::seq(lo, hi, N, endpoints);
	return Rcpp::NumericVector(out.begin(), out.end());
}

