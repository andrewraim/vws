#ifndef VWS_CATEG_H
#define VWS_CATEG_H

#include <Rcpp.h>
#include "gumbel.h"

namespace vws {

inline unsigned int r_categ(const Rcpp::NumericVector& p, bool log = false,
	bool one_based = false)
{
	unsigned int k = p.length();
	Rcpp::NumericVector lp;
	if (log) {
		lp = p;
	} else {
		lp = Rcpp::log(p);
	}

	const Rcpp::NumericVector& z = r_gumbel(k);
	return Rcpp::which_max(z + lp) + one_based;
}

inline Rcpp::IntegerVector r_categ(unsigned int n,
	const Rcpp::NumericVector& p, bool log = false, bool one_based = false)
{
	Rcpp::IntegerVector out(n);

	for (unsigned int i = 0; i < n; i++) {
		out(i) = r_categ(p, log, one_based);
	}

	return out;
}

}

#endif
