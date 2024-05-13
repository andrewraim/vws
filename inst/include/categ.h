#ifndef VWS_CATEG_H
#define VWS_CATEG_H

#include <Rcpp.h>
#include "gumbel.h"

namespace vws {

unsigned int r_categ(const Rcpp::NumericVector& p, bool log = false)
{
	unsigned int k = p.length();
	Rcpp::NumericVector lp;
	if (log) {
		lp = p;
	} else {
		lp = Rcpp::log(p);
	}

	const Rcpp::NumericVector& z = r_gumbel(k);
	return Rcpp::which_max(z + lp);		// TBD: are these zero-based indices?
}

Rcpp::IntegerVector r_categ(unsigned int n, const Rcpp::NumericVector& p, bool log = false)
{
	Rcpp::IntegerVector out(n);

	for (unsigned int i = 0; i < n; i++) {
		out(i) = r_categ(p, log);
	}

	return out;
}

}

#endif
