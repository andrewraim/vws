#ifndef RCPP_CATEG_H
#define RCPP_CATEG_H

#include <Rcpp.h>

//' Categorical Distribution
//'
//' @name Categorical
//' @export
// [[Rcpp::export(name = "r_categ")]]
Rcpp::IntegerVector r_categ_rcpp(unsigned int n,
	const Rcpp::NumericVector& p, bool log = false, bool one_based = false);

#endif

