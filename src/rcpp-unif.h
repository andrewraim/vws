#ifndef RCPP_UNIF_H
#define RCPP_UNIF_H

#include <Rcpp.h>

//' Uniform Distribution
//'
//' @name Uniform
//' @export
// [[Rcpp::export(name = "mgf_unif")]]
Rcpp::NumericVector mgf_unif_rcpp(const Rcpp::NumericVector& s, double lo,
	double hi, bool log = false);

#endif
