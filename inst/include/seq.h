#ifndef SEQ_H
#define SEQ_H

#include <Rcpp.h>

/*
* Function to generate a sequence of numbers.
*/

namespace vws {

/*
* Prepare knots which define $N$ equally-spaced intervals between (finite)
* endpoints `lo` and `hi`.
*
* - `lo`: Left endpoint; must be finite.
* - `hi`: Right endpoint; must be finite.
* - `N`: Number of desired intervals.
* - `endpoints`: logical; if `true`, include endpoints
*/
inline std::vector<double> seq(double lo, double hi, unsigned int N, bool endpoints = false)
{
	if (N == 0) { Rcpp::stop("N must be positive"); }
	if (std::isinf(lo)) { Rcpp::stop("lo must be finite"); };
	if (std::isinf(hi)) { Rcpp::stop("hi must be finite"); };
	if (hi < lo) { Rcpp::stop("Require lo < hi"); };


	std::vector<double> out;
	if (endpoints) { out.push_back(lo); }

	double x = lo;
	double delta = (hi - lo) / N;

	for (unsigned int i = 0; i < N-1; i++) {
		x += delta;
		out.push_back(x);
	}

	if (endpoints) { out.push_back(hi); }

	return out;
}

}

#endif
