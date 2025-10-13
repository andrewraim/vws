#ifndef VWS_TYPEDEFS_H
#define VWS_TYPEDEFS_H

#include <Rcpp.h>
#include "fntl.h"
#include "optimize-hybrid.h"

namespace vws {

/*
* A constant that represents the maximum value of unsigned int.
*/
static const unsigned int uint_max = std::numeric_limits<unsigned int>::max();

/*
* Typedefs for several functions used in VWS programming.
*/
typedef std::function<double(double, bool)> dfdb;
typedef std::function<bool(double)> indicator;
typedef std::function<double(const dfdb& w, double lo, double hi,
	bool log)> optimizer;

/*
* Enumeration that represents actions to be taken when an error condition
* arises. The general meanings are:
*
* - `STOP`: an exception should be thrown;
* - `WARNING`: a warning should be emitted but execution should continue;
* - `MESSAGE`: a message should be emitted but execution should continue;
* - `NONE`: the condition should be ignored.
*/
enum class error_action {
	STOP,
	WARNING,
	MESSAGE,
	NONE
};

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
