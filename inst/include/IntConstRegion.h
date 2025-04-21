#ifndef VWS_INT_CONST_REGION_H
#define VWS_INT_CONST_REGION_H

#include <Rcpp.h>
#include <memory>
#include "fntl.h"
#include "Region.h"
#include "RealConstRegion.h"
#include "UnivariateHelper.h"
#include "optimize-hybrid.h"
#include "log-sum-exp.h"
#include "typedefs.h"

namespace vws {

/*
* Integer-Valued Region with Constant Majorizer
*
* A subclass of Region based on univariate intervals and a constant majorizer
* for the weight function. This version is for integer supports. It uses
* doubles for its domain but values are integers.
*/
class IntConstRegion : public RealConstRegion
{
public:
	/*
	* Construct a singleton region based on interval $(a,a]$.
	* - `a`: Lower and upper limit of interval.
	* - `w`: Weight function for the target distribution.
	* - `helper`: contains operations of the base distribution $g$.
	*/
	IntConstRegion(double a, const uv_weight_function& w,
		const UnivariateHelper& helper);

	/*
	* Construct a region based on interval $(a,b]$.
	* - `a` Lower limit of interval.
	* - `b` Upper limit of interval.
	* - `w` Weight function for the target distribution.
	* - `helper`: contains operations of the base distribution $g$.
	*/
	IntConstRegion(double a, double b, const uv_weight_function& w,
		const UnivariateHelper& helper);

	/*
	* The following functions override methods in `RealConstRegion`. See that
	* class' documentation for their interfaces.
	*/
	// bool s(const double& x) const;
	bool is_bifurcatable() const;

	/*
	* A midpoint between limits $a$ and $b$ of region. If $a$ and $b$ are both
	* finite, return the standard midpoint. If both are infinite, return zero.
	* If only $a$ is finite, return a larger point in the support. If only $b$
	* is finite, return a smaller point in the support.
	*/
	// double midpoint() const;

	/*
	* Return a pair of regions that result from bifurcating this region. The
	* bifurcation point is chosen to be the midpoint of $(a, b]$.
	*/
	std::pair<IntConstRegion,IntConstRegion> bifurcate() const;

	/*
	* Return a pair of regions that result from bifurcating this region at $x$.
	*/
	std::pair<IntConstRegion,IntConstRegion> bifurcate(const double& x) const;

	/*
	* Return a region based on the singleton interval $(x, x]$, using this
	* object's weight function, base distribution, etc.
	*/
	IntConstRegion singleton(const double& x) const;

	/*
	* See RealConstRegion for documentation of the following.
	*/
	bool operator<(const IntConstRegion& x) const;
	bool operator==(const IntConstRegion& x) const;
	const IntConstRegion& operator=(const IntConstRegion& x);
};

inline IntConstRegion::IntConstRegion(double a,
	const uv_weight_function& w, const UnivariateHelper& helper)
: RealConstRegion(a, w, helper)
{
}

inline IntConstRegion::IntConstRegion(double a, double b,
	const uv_weight_function& w, const UnivariateHelper& helper)
: RealConstRegion(a, b, w, helper)
{
}

/*
inline bool IntConstRegion::s(const double& x) const
{
	return (_a < x && x <= _b) && _helper->s(x);
}
*/

/*
inline double IntConstRegion::midpoint() const
{
	Rcpp::stop("Implement me");

	double out;

	if (std::isinf(_a) && std::isinf(_b) && _a < 0 && _b > 0) {
		// In this case, we have an interval (-inf, inf). Make a split at zero.
		out = 0;
	} else if (std::isinf(_a) && _a < 0) {
		// Left endpoint is -inf. Split based on right endpoint.
		out = _b - std::fabs(_b) - 1;
	} else if (std::isinf(_b) && _b > 0) {
		// Right endpoint is inf. Split based on left endpoint.
		out = _a + std::fabs(_a) + 1;
	} else {
		out = (_a + _b) / 2;
	}

	return out;
}
*/

inline std::pair<IntConstRegion,IntConstRegion>
IntConstRegion::bifurcate() const
{
	return bifurcate(midpoint());
}

inline std::pair<IntConstRegion,IntConstRegion>
IntConstRegion::bifurcate(const double& x) const
{
	IntConstRegion r1(_a, x, *_w, *_helper);
	IntConstRegion r2(x, _b, *_w, *_helper);
	return std::make_pair(r1, r2);
}

inline IntConstRegion IntConstRegion::singleton(const double& x) const
{
	return IntConstRegion(x, *_w, *_helper);
}

inline bool IntConstRegion::is_bifurcatable() const
{
	// Return true if there are at least two integers between a and b
	// Rprintf("Region (%g,%g] contains %g integers. is_bifurcatable = %d\n",
	// 	_a, _b, std::floor(_b) - std::ceil(_a) + 1,
	// 	std::ceil(_a) + 1 <= std::floor(_b));
	// return std::ceil(_a) + 1 <= std::floor(_b);
	return _b - _a > 1;
}

inline bool IntConstRegion::operator<(const IntConstRegion& x) const
{
	// return _a < x._a;
	return _b <= x._a;
}

inline bool IntConstRegion::operator==(const IntConstRegion& x) const
{
	return _a == x._a && _b == x._b;
}

inline const IntConstRegion& IntConstRegion::operator=(const IntConstRegion& x)
{
	_a = x._a;
	_b = x._b;
	_w = x._w;
	_helper = x._helper;
	_log_w_max = x._log_w_max;
	_log_w_min = x._log_w_min;
	_log_prob = x._log_prob;
	return *this;
}

}

#endif
