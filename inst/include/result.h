#ifndef VWS_RESULT_H
#define VWS_RESULT_H

#include <RcppCommon.h>
#include <vector>

namespace vws {

template <typename T>
struct rejection_result
{
	std::vector<T> draws;
	std::vector<unsigned int> rejects;

	operator SEXP() const;
};

}

#include <Rcpp.h>

namespace vws {

template <typename T>
inline rejection_result<T>::operator SEXP() const
{
	return Rcpp::List::create(
		Rcpp::Named("draws") = draws,
		Rcpp::Named("rejects") = rejects
	);
}

}

#endif
