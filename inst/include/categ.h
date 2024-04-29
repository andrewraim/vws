#ifndef CATEG_H
#define CATEG_H

#include <Rcpp.h>

double r_categ(unsigned int n, const Rcpp::NumericVector& p, bool log_p = false)
{
	unsigned int k = length(p);
	double lp = log_p ? p : log(p);
	double z = r_gumbel(k);
	return Rcpp::which_max(z + lp);		// TBD: are these zero-based indices?
}

#endif
