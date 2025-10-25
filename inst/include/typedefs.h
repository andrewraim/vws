#ifndef VWS_TYPEDEFS_H
#define VWS_TYPEDEFS_H

#include <Rcpp.h>
#include "fntl.h"
#include "optimize-hybrid.h"

namespace vws {

/*
* Typedefs for several functions used in VWS programming.
*/
typedef std::function<double(double, bool)> dfdb;
typedef std::function<double(const dfdb& w, double lo, double hi,
	bool log)> optimizer;

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

}

#endif

