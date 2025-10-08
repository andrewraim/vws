#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <Rcpp.h>

// MGF needed for choosing expansion point in LinearVWSRegion
// [[Rcpp::export]]
double mgf_truncpois(double s, double lambda, double a, double b,
	bool log = false);

// [[Rcpp::export]]
Rcpp::NumericVector incgamma(double a, const Rcpp::NumericVector& x,
	bool lower = false, bool log = false);

double incgamma(double a, double x, bool lower = false, bool log = false);

#endif
