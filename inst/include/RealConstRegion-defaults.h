#ifndef VWS_REAL_CONST_REGION_DEFAULTS_H
#define VWS_REAL_CONST_REGION_DEFAULTS_H

#include <Rcpp.h>
#include "typedefs.h"

namespace vws {

inline static const optimizer maxopt_default = [](const dfdb& w,
	double lo, double hi, bool log) -> double
{
	// Pass the log-weight function to `optimize_hybrid`.
	const fntl::dfd& f = [&](double x) -> double { return w(x, true); };
	const auto& out = optimize_hybrid(f, 0, lo, hi, true);
	return log ? out.value : exp(out.value);
};

inline static const optimizer minopt_default = [](const dfdb& w,
	double lo, double hi, bool log) -> double
{
	// Pass the log-weight function to `optimize_hybrid`.
	const fntl::dfd& f = [&](double x) -> double { return w(x, true); };
	const auto& out = optimize_hybrid(f, 0, lo, hi, false);
	return log ? out.value : exp(out.value);
};

/*
* A modified version of the arithmetic midpoint. When one endpoint is infinite,
* bring in the finite endpoint by a factor of 1/2 and an additional unit of one.
* This helps to avoid large numbers of wasted regions during knot selection.
* When both endpoints are infinite, take zero to be the midpoint.
*/
inline static const midpoint midpoint_default = [](double a, double b) -> double
{
	double out;

	if (std::isinf(a) && std::isinf(b) && a < 0 && b > 0) {
		// In this case, we have an interval (-inf, inf). Make a split at zero.
		out = 0;
	} else if (std::isinf(a) && a < 0) {
		// Left endpoint is -inf. Split based on right endpoint.
		double sgn = (b > 0) - (b < 0);
		out = b * std::pow(2, -sgn) - 1;
	} else if (std::isinf(b) && b > 0) {
		// Right endpoint is inf. Split based on left endpoint.
		double sgn = (a > 0) - (a < 0);
		out = a * std::pow(2, sgn) + 1;
	} else {
		out = (a + b) / 2;
	}

	return out;
};

}

#endif

