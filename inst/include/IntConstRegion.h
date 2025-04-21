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
	* - `opt`: a function of type `weight_optimization` that can both minimize
	*   and maximize `w`. If not specified, we use numerical optimization.
	*/
	IntConstRegion(double a, const uv_weight_function& w,
		const UnivariateHelper& helper);
	IntConstRegion(double a, const uv_weight_function& w,
		const UnivariateHelper& helper, const weight_optimization& opt);

	/*
	* Construct a region based on interval $(a,b]$.
	* - `a` Lower limit of interval.
	* - `b` Upper limit of interval.
	* - `w` Weight function for the target distribution.
	* - `helper`: contains operations of the base distribution $g$.
	* - `opt`: a function of type `weight_optimization` that can both minimize
	*   and maximize `w`. If not specified, we use numerical optimization.
	*/
	IntConstRegion(double a, double b, const uv_weight_function& w,
		const UnivariateHelper& helper);
	IntConstRegion(double a, double b, const uv_weight_function& w,
		const UnivariateHelper& helper, const weight_optimization& opt);

	/*
	* The following functions override methods in `RealConstRegion`. See that
	* class' documentation for their interfaces.
	*/
	// bool s(const double& x) const;
	bool is_bifurcatable() const;

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

inline IntConstRegion::IntConstRegion(double a,
	const uv_weight_function& w, const UnivariateHelper& helper,
	const weight_optimization& opt)
: RealConstRegion(a, w, helper, opt)
{
}

inline IntConstRegion::IntConstRegion(double a, double b,
	const uv_weight_function& w, const UnivariateHelper& helper)
: RealConstRegion(a, b, w, helper)
{
}

inline IntConstRegion::IntConstRegion(double a, double b,
	const uv_weight_function& w, const UnivariateHelper& helper,
	const weight_optimization& opt)
: RealConstRegion(a, b, w, helper, opt)
{
}

inline std::pair<IntConstRegion,IntConstRegion>
IntConstRegion::bifurcate() const
{
	return bifurcate(midpoint());
}

inline std::pair<IntConstRegion,IntConstRegion>
IntConstRegion::bifurcate(const double& x) const
{
	IntConstRegion r1(_a, x, *_w, *_helper, _opt);
	IntConstRegion r2(x, _b, *_w, *_helper, _opt);
	return std::make_pair(r1, r2);
}

inline IntConstRegion IntConstRegion::singleton(const double& x) const
{
	return IntConstRegion(x, *_w, *_helper, _opt);
}

inline bool IntConstRegion::is_bifurcatable() const
{
	// Return true if the distance between a and b allows for two or more
	// integers
	return _b - _a >= 2;
}

inline bool IntConstRegion::operator<(const IntConstRegion& x) const
{
	return RealConstRegion::operator<(x);
}

inline bool IntConstRegion::operator==(const IntConstRegion& x) const
{
	return RealConstRegion::operator==(x);
}

inline const IntConstRegion& IntConstRegion::operator=(const IntConstRegion& x)
{
	RealConstRegion::operator=(x);
	return *this;
}

}

#endif
