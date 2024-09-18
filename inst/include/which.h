#ifndef VWS_WHICH_H
#define VWS_WHICH_H

#include <Rcpp.h>

namespace vws {

/*
* A C++ implementation of the R "which" function. Does a standard
* implementation already exist?
*
* Returns the indices in `x` which have the value `true`. If `one_based` is
* true, indices are transformed to start at 1 (especially for use in R).
* Otherwise, indices begin at 0.
*/
inline Rcpp::IntegerVector which(const Rcpp::LogicalVector& x, bool one_based = false)
{
	std::vector<unsigned int> idx;

	for (unsigned int i = 0; i < x.size(); i++) {
		if (x(i)) {
			idx.push_back(i + one_based);
		}
	}

	Rcpp::IntegerVector out(idx.begin(), idx.end());
	return out;
}

}

#endif
