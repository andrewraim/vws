#ifndef VWS_UNIF_H
#define VWS_UNIF_H

#include <Rcpp.h>
#include "log-sum-exp.h"

namespace vws {

inline double mgf_unif(double s, double lo, double hi, bool log)
{
	double out;

	if (s < 0) {
		out = log_sub2_exp(lo*s, hi*s) - std::log(hi - lo) - std::log(-s);
	} else if (s > 0) {
		out = log_sub2_exp(hi*s, lo*s) - std::log(hi - lo) - std::log(s);
	} else {
		out = R_NegInf;
	}

	return log ? out : std::exp(out);
}

inline Rcpp::NumericVector mgf_unif(const Rcpp::NumericVector& s, double lo,
	double hi, bool log = false)
{
	unsigned int n = s.length();
	Rcpp::NumericVector out(n);

	for (unsigned int i = 0; i < n; i++) {
		out(i) = mgf_unif(s(i), lo, hi, log);
	}

	return out;
}

}

#endif
