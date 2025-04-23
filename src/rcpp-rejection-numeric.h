#ifndef RCPP_REJECTION_NUMERIC_H
#define RCPP_REJECTION_NUMERIC_H

#include <Rcpp.h>

// [[Rcpp::export]]
Rcpp::List rejection_numeric_rcpp(unsigned int n, double lo, double hi,
	const Rcpp::Function& w, const Rcpp::Function& d_base,
	const Rcpp::Function& p_base, const Rcpp::Function& q_base,
	unsigned int N, double tol, const Rcpp::List& control);

// [[Rcpp::export]]
Rcpp::List rejection_numeric_opt_rcpp(unsigned int n, double lo, double hi,
	const Rcpp::Function& w, const Rcpp::Function& d_base,
	const Rcpp::Function& p_base, const Rcpp::Function& q_base,
	const Rcpp::Function& maxopt, const Rcpp::Function& minopt,
	unsigned int N, double tol, const Rcpp::List& control);

#endif

