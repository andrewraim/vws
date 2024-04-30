#ifndef CATEG_H
#define CATEG_H

#include <Rcpp.h>
#include "gumbel.h"

namespace vws {

double r_categ(unsigned int n, const Rcpp::NumericVector& p, bool log = false)
{
	unsigned int k = p.length();
	// const Rcpp::NumericVector& lp = log ? p : Rcpp::log(p);
	Rcpp::NumericVector lp;
	if (log) {
		lp = p;
	} else {
		lp = Rcpp::log(p);
	}

	const Rcpp::NumericVector& z = r_gumbel(k);
	return Rcpp::which_max(z + lp);		// TBD: are these zero-based indices?
}

}

#endif
