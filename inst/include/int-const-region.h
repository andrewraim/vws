#ifndef VWS_INT_CONST_REGION_H
#define VWS_INT_CONST_REGION_H

#include <Rcpp.h>
#include "typedefs.h"
#include "real-const-region.h"
#include "real-const-region-defaults.h"
#include "univariate-helper.h"

namespace vws {

/*
* Integer-Valued Region with Constant Majorizer
*
* A subclass of Region based on univariate intervals and a constant majorizer
* for the weight function. This version is for integer supports. It uses
* doubles for its domain but values are integers.
*/
class int_const_region : public real_const_region
{
public:
	/*
	* Construct a singleton region based on interval $(a,a]$.
	* - `a`: Lower and upper limit of interval.
	* - `w`: Weight function for the target distribution.
	* - `helper`: contains operations of the base distribution $g$.
	* - `maxopt`: a function of type `optimizer` that maximizes `w`.
	* - `minopt`: a function of type `optimizer` that minimizes `w`.
	* - `mid`: a function of type `midpoint` to compute the midpoint of interval
	*   regions.
	*
	* If `maxopt` and `minopt` are not specified, we use numerical optimization.
	* If `mid` is not specified, we use a version of the arithmetic midpoint
	* with special handling of infinite endpoints.
	*/
	int_const_region(double a, const dfdb& w,
		const univariate_helper& helper,
		const optimizer& maxopt = maxopt_default,
		const optimizer& minopt = minopt_default,
		const vws::midpoint& mid = midpoint_default);

	/*
	* Construct a region based on interval $(a,b]$.
	* - `a` Lower limit of interval.
	* - `b` Upper limit of interval.
	* - `w` Weight function for the target distribution.
	* - `helper`: contains operations of the base distribution $g$.
	* - `maxopt`: a function of type `optimizer` that maximizes `w`.
	* - `minopt`: a function of type `optimizer` that minimizes `w`.
	* - `mid`: a function of type `midpoint` to compute the midpoint of interval
	*   regions.
	*
	* If `maxopt` and `minopt` are not specified, we use numerical optimization.
	* If `mid` is not specified, we use a version of the arithmetic midpoint
	* with special handling of infinite endpoints.
	*/
	int_const_region(double a, double b, const dfdb& w,
		const univariate_helper& helper,
		const optimizer& maxopt = maxopt_default,
		const optimizer& minopt = minopt_default,
		const vws::midpoint& mid = midpoint_default);

	/*
	* The following functions override methods in `real_const_region`. See that
	* class' documentation for their interfaces.
	*/
	bool is_bifurcatable() const;

	/*
	* Return a pair of regions that result from bifurcating this region. The
	* bifurcation point is chosen to be the midpoint of $(a, b]$.
	*/
	std::pair<int_const_region,int_const_region> bifurcate() const;

	/*
	* Return a pair of regions that result from bifurcating this region at $x$.
	*/
	std::pair<int_const_region,int_const_region> bifurcate(const double& x) const;

	/*
	* Return a region based on the singleton interval $(x, x]$, using this
	* object's weight function, base distribution, etc.
	*/
	int_const_region singleton(const double& x) const;

	/*
	* See real_const_region for documentation of the following.
	*/
	bool operator<(const int_const_region& x) const;
	bool operator==(const int_const_region& x) const;
	const int_const_region& operator=(const int_const_region& x);
};

inline int_const_region::int_const_region(double a,
	const dfdb& w, const univariate_helper& helper,
	const optimizer& maxopt, const optimizer& minopt, const vws::midpoint& mid)
: real_const_region(a, w, helper, maxopt, minopt, mid)
{
}

inline int_const_region::int_const_region(double a, double b,
	const dfdb& w, const univariate_helper& helper,
	const optimizer& maxopt, const optimizer& minopt, const vws::midpoint& mid)
: real_const_region(a, b, w, helper, maxopt, minopt, mid)
{
}

inline std::pair<int_const_region,int_const_region>
int_const_region::bifurcate() const
{
	return bifurcate(midpoint());
}

inline std::pair<int_const_region,int_const_region>
int_const_region::bifurcate(const double& x) const
{
	int_const_region r1(_a, x, _w, _helper, _maxopt, _minopt, _mid);
	int_const_region r2(x, _b, _w, _helper, _maxopt, _minopt, _mid);
	return std::make_pair(r1, r2);
}

inline int_const_region int_const_region::singleton(const double& x) const
{
	return int_const_region(x, _w, _helper, _maxopt, _minopt, _mid);
}

inline bool int_const_region::is_bifurcatable() const
{
	// Return true if the distance between a and b allows for two or more
	// integers
	return _b - _a > 1;
}

inline bool int_const_region::operator<(const int_const_region& x) const
{
	return real_const_region::operator<(x);
}

inline bool int_const_region::operator==(const int_const_region& x) const
{
	return real_const_region::operator==(x);
}

inline const int_const_region& int_const_region::operator=(const int_const_region& x)
{
	real_const_region::operator=(x);
	return *this;
}

}

#endif

