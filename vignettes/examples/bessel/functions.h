#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <Rcpp.h>

// MGF needed for choosing expansion point in LinearVWSRegion
// [[Rcpp::export]]
double mgf_truncpois(double s, double a, double b, double lambda,
	bool log = false);

// [[Rcpp::export]]
double incgamma(double a, double x, bool lower = false, bool log = false);

#endif
