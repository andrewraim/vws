#ifndef VWS_LOGIT_H
#define VWS_LOGIT_H

#include <Rcpp.h>

namespace vws {

inline double logit(double p)
{
	return R::plogis(p, 0, 1, true, false);
}

inline double inv_logit(double x)
{
	return R::qlogis(x, 0, 1, true, false);
}

inline Rcpp::NumericVector logit(const Rcpp::NumericVector& p)
{
	return Rcpp::plogis(p);
}

inline Rcpp::NumericVector inv_logit(const Rcpp::NumericVector& x)
{
	return Rcpp::qlogis(x);
}

}

#endif
