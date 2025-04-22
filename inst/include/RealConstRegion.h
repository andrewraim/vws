#ifndef VWS_REAL_CONST_REGION_H
#define VWS_REAL_CONST_REGION_H

#include <Rcpp.h>
#include <memory>
#include "fntl.h"
#include "Region.h"
#include "UnivariateHelper.h"
#include "optimize-hybrid.h"
#include "log-sum-exp.h"
#include "typedefs.h"

namespace vws {

/*
* Real-Valued Region with Constant Majorizer
*
* A subclass of Region based on univariate intervals and a constant majorizer
* for the weight function. This version is for continuous supports.
*/
class RealConstRegion : public Region<double>
{
public:
	/*
	* Construct a singleton region based on interval $(a,a]$.
	* - `a`: Lower and upper limit of interval.
	* - `w`: Weight function for the target distribution.
	* - `helper`: contains operations of the base distribution $g$.
    * - `maxopt`: a function of type `optimizer` that maximizes `w`.
    * - `minopt`: a function of type `optimizer` that minimizes `w`.
    *
    * If `maxopt` and `minopt` are not specified, we use numerical optimization.
    *
	*/
	RealConstRegion(double a, const uv_weight_function& w,
		const UnivariateHelper& helper);
	RealConstRegion(double a, const uv_weight_function& w,
		const UnivariateHelper& helper, const optimizer& maxopt,
		const optimizer& minopt);

	/*
	* Construct a region based on interval $(a,b]$.
	* - `a` Lower limit of interval.
	* - `b` Upper limit of interval.
	* - `w` Weight function for the target distribution.
	* - `helper`: contains operations of the base distribution $g$.
    * - `maxopt`: a function of type `optimizer` that maximizes `w`.
    * - `minopt`: a function of type `optimizer` that minimizes `w`.
    *
    * If `maxopt` and `minopt` are not specified, we use numerical optimization.
	*/
	RealConstRegion(double a, double b, const uv_weight_function& w,
		const UnivariateHelper& helper);
	RealConstRegion(double a, double b, const uv_weight_function& w,
		const UnivariateHelper& helper, const optimizer& maxopt,
		const optimizer& minopt);

	/*
	* The following functions override abstract methods in `Region`. See that
	* class' documentation for their interfaces.
	*/
	double d_base(const double& x, bool log = false) const;
	std::vector<double> r(unsigned int n) const;
	bool s(const double& x) const;
	double w(const double& x, bool log = true) const;
	double w_major(const double& x, bool log = true) const;
	bool is_bifurcatable() const;
	double get_lower() const { return _a; }
	double get_upper() const { return _b; }
	double get_xi_upper(bool log = true) const;
	double get_xi_lower(bool log = true) const;
	std::string description() const;

	/*
	* Maximize or minimize the function $w$ over this region. Optimization
	* is carried out with the `optimize_hybrid` function.
	*
	* - `maximize`: if `true` do maximization; otherwise do minimization.
	* - `log`: if `true`, return value on the log-scale. Otherwise, return it
	*   on the original scale.
	*
	* Returns the optimized value of $w$.
	*/
	double optimize(bool maximize = true, bool log = true) const;

	/*
	* A midpoint between limits $a$ and $b$ of region. If $a$ and $b$ are both
	* finite, return the standard midpoint. If both are infinite, return zero.
	* If only $a$ is finite, return a larger point in the support. If only $b$
	* is finite, return a smaller point in the support.
	*/
	double midpoint() const;

	/*
	* Return a pair of regions that result from bifurcating this region. The
	* bifurcation point is chosen to be the midpoint of $(a, b]$.
	*/
	std::pair<RealConstRegion,RealConstRegion> bifurcate() const;

	/*
	* Return a pair of regions that result from bifurcating this region at $x$.
	*/
	std::pair<RealConstRegion,RealConstRegion> bifurcate(const double& x) const;

	/*
	* Return a region based on the singleton interval $(x, x]$, using this
	* object's weight function, base distribution, etc.
	*/
	RealConstRegion singleton(const double& x) const;

	/*
	* Region $(a_1, b_1]$ is considered "less than" $(a_2, b_2]$ if $a_1 < a_2$.
	*/
	bool operator<(const RealConstRegion& x) const;

	/*
	* Region $(a_1, b_1]$ is considered "equal to" $(a_2, b_2]$ if $a_1 = a_2$
	* and $b_1 = b_2$.
	*/
	bool operator==(const RealConstRegion& x) const;

	/*
	* Set this Region to be equal to `x`.
	*/
	const RealConstRegion& operator=(const RealConstRegion& x);

protected:
	void init();
	double _a;
	double _b;
	const uv_weight_function* _w;
	const UnivariateHelper* _helper;
	double _log_w_max;
	double _log_w_min;
	double _log_prob;
	optimizer _maxopt;
	optimizer _minopt;
};

inline RealConstRegion::RealConstRegion(double a,
	const uv_weight_function& w, const UnivariateHelper& helper)
: _a(a), _b(a), _w(&w), _helper(&helper), _log_w_max(NAN), _log_w_min(NAN),
  _log_prob(NAN), _maxopt(), _minopt()
{
	_log_w_max = (*_w)(a, true);
	_log_w_min = (*_w)(a, true);
	_log_prob = R_NegInf;
}

inline RealConstRegion::RealConstRegion(double a,
	const uv_weight_function& w, const UnivariateHelper& helper,
	const optimizer& maxopt, const optimizer& minopt)
: _a(a), _b(a), _w(&w), _helper(&helper), _log_w_max(NAN), _log_w_min(NAN),
  _log_prob(NAN), _maxopt(maxopt), _minopt(minopt)
{
	RealConstRegion::init();
}

inline RealConstRegion::RealConstRegion(double a, double b,
	const uv_weight_function& w, const UnivariateHelper& helper)
: _a(a), _b(b), _w(&w), _helper(&helper), _log_w_max(NAN), _log_w_min(NAN),
  _log_prob(NAN), _maxopt(), _minopt()
{
	RealConstRegion::init();
}

inline RealConstRegion::RealConstRegion(double a, double b,
	const uv_weight_function& w, const UnivariateHelper& helper,
	const optimizer& maxopt, const optimizer& minopt)
: _a(a), _b(b), _w(&w), _helper(&helper), _log_w_max(NAN), _log_w_min(NAN),
  _log_prob(NAN), _maxopt(maxopt), _minopt(minopt)
{
	RealConstRegion::init();
}

inline void RealConstRegion::init()
{
	if (_a > _b) {
		// Invalid interval
		Rcpp::stop("a > b");
	} else if (_a < _b) {
		if (_maxopt) {
			// Use custom optimization function that was provided
			_log_w_max = _maxopt(*_w, _a, _b, true);
		} else {
			// Use default numerical optimization
			_log_w_max = optimize(true);
		}

		if (_minopt) {
			// Use custom optimization function that was provided
			_log_w_min = _minopt(*_w, _a, _b, true);
		} else {
			// Use default numerical optimization
			_log_w_min = optimize(false);
		}

		// Compute $P(a < X <= b)$ for $X \sim g$ on the log scale.
		_log_prob = log_sub2_exp(_helper->p(_b, true, true), _helper->p(_a, true, true));
	} else {
		// Singleton interval (a,a]
		_log_w_max = (*_w)(_a, true);
		_log_w_min = (*_w)(_a, true);
		_log_prob = R_NegInf;
	}

	if ( (_log_w_max > 0) && std::isinf(_log_w_max) ) {
		Rcpp::stop("%s. %s",
			"Infinite maximum value found in optimize",
			"Cannot be used with RealConstRegion");
	}
}

inline double RealConstRegion::d_base(const double& x, bool log) const
{
	return _helper->d(x, log);
}

inline std::vector<double> RealConstRegion::r(unsigned int n) const
{
	// Generate a draw from $g_j$; i.e., the density $g$ truncated to this
	// region. Compute g$q((pb - pa) * u + pa) on the log scale.
	const Rcpp::NumericVector& u = Rcpp::runif(n);
	double log_pa = _helper->p(_a, true, true);
	const Rcpp::NumericVector& log_p = log_add2_exp(_log_prob + log(u), Rcpp::rep(log_pa, n));

	std::vector<double> out;
	for (unsigned int i = 0; i < n; i++) {
		out.push_back(_helper->q(log_p(i), true, true));
	}

	return out;
}

inline bool RealConstRegion::s(const double& x) const
{
	return (_a < x && x <= _b) && _helper->s(x);
}

inline double RealConstRegion::w(const double& x, bool log) const
{
	return (*_w)(x, log);
}

inline double RealConstRegion::w_major(const double& x, bool log) const
{
	double out = _helper->s(x) ? _log_w_max : R_NegInf;
	return log ? out : exp(out);
}

inline double RealConstRegion::midpoint() const
{
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

inline std::pair<RealConstRegion,RealConstRegion>
RealConstRegion::bifurcate() const
{
	return bifurcate(midpoint());
}

inline std::pair<RealConstRegion,RealConstRegion>
RealConstRegion::bifurcate(const double& x) const
{
	RealConstRegion r1(_a, x, *_w, *_helper, _maxopt, _minopt);
	RealConstRegion r2(x, _b, *_w, *_helper, _maxopt, _minopt);
	return std::make_pair(r1, r2);
}

inline RealConstRegion RealConstRegion::singleton(const double& x) const
{
	return RealConstRegion(x, *_w, *_helper, _maxopt, _minopt);
}

inline bool RealConstRegion::is_bifurcatable() const
{
	return true;
}

inline double RealConstRegion::get_xi_upper(bool log) const
{
	double log_xi_upper = _log_w_max + _log_prob;
	return log ? log_xi_upper : exp(log_xi_upper);
}

inline double RealConstRegion::get_xi_lower(bool log) const
{
	double log_xi_lower = _log_w_min + _log_prob;
	return log ? log_xi_lower : exp(log_xi_lower);
}

inline std::string RealConstRegion::description() const
{
	char buf[32];
	sprintf(buf, "(%g, %g]", _a, _b);
	return buf;
}

inline double RealConstRegion::optimize(bool maximize, bool log) const
{
	// Pass the log-weight function to `optimize_hybrid`.
    const fntl::dfd& f = [&](double x) -> double { return (*_w)(x, true); };
	const auto& out = optimize_hybrid(f, 0, _a, _b, maximize);
	return log ? out.value : exp(out.value);
}

inline bool RealConstRegion::operator<(const RealConstRegion& x) const
{
	return _b <= x._a;
}

inline bool RealConstRegion::operator==(const RealConstRegion& x) const
{
	return _a == x._a && _b == x._b;
}

inline const RealConstRegion& RealConstRegion::operator=(const RealConstRegion& x)
{
	_a = x._a;
	_b = x._b;
	_w = x._w;
	_helper = x._helper;
	_log_w_max = x._log_w_max;
	_log_w_min = x._log_w_min;
	_log_prob = x._log_prob;
	_maxopt = x._maxopt;
	_minopt = x._minopt;
	return *this;
}

}

#endif
