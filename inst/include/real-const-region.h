#ifndef VWS_REAL_CONST_REGION_H
#define VWS_REAL_CONST_REGION_H

#include <Rcpp.h>
#include "typedefs.h"
#include "region.h"
#include "real-const-region-defaults.h"
#include "univariate-helper.h"
#include "log-sum-exp.h"

namespace vws {

/*
* Real-Valued region with Constant Majorizer
*
* A subclass of region based on univariate intervals and a constant majorizer
* for the weight function. This version is for continuous supports.
*/
class real_const_region : public region<double>
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
	real_const_region(double a, const dfdb& w,
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
	*
	* If `maxopt` and `minopt` are not specified, we use numerical optimization.
	* If `mid` is not specified, we use a version of the arithmetic midpoint
	* with special handling of infinite endpoints.
	*/
	real_const_region(double a, double b, const dfdb& w,
		const univariate_helper& helper,
		const optimizer& maxopt = maxopt_default,
		const optimizer& minopt = minopt_default,
		const vws::midpoint& mid = midpoint_default);

	/*
	* The following functions override abstract methods in `region`. See that
	* class' documentation for their interfaces.
	*/
	double d_base(const double& x, bool log = false) const;
	std::vector<double> r(unsigned int n) const;
	bool s(const double& x) const;
	double w(const double& x, bool log = true) const;
	double w_major(const double& x, bool log = true) const;
	double w_minor(const double& x, bool log = true) const;
	bool is_bifurcatable() const;
	bool is_mergeable(const real_const_region& x) const;
	double lower() const { return _a; }
	double upper() const { return _b; }
	double xi_upper(bool log = true) const;
	double xi_lower(bool log = true) const;
	std::string description() const;

	/*
	* Set methods for updating existing regions.
	*
	* Important: the init function should be called after using set methods and
	* before the region is used for sampling. The init method recomputes saved
	* quantities such as the majorization and minorization constants needed for
	* rejection sampling.
	*/
	void set_w(const dfdb& w);
	void set_helper(const univariate_helper& helper);
	void set_maxopt(const optimizer& maxopt);
	void set_minopt(const optimizer& minopt);
	void set_mid(const vws::midpoint& mid);
	void init();

	/*
	* Return a region that consists of this one merged with argument `x`.
	*/
	real_const_region merge(const real_const_region& x) const;

	/*
	* Return the midpoint of this region using the function _mid.
	*/
	double midpoint() const;

	/*
	* Return a pair of regions that result from bifurcating this region. The
	* bifurcation point is chosen to be the midpoint of $(a, b]$.
	*/
	std::pair<real_const_region,real_const_region> bifurcate() const;

	/*
	* Return a pair of regions that result from bifurcating this region at $x$.
	*/
	std::pair<real_const_region,real_const_region> bifurcate(const double& x) const;

	/*
	* Return a region based on the singleton interval $(x, x]$, using this
	* object's weight function, base distribution, etc.
	*/
	real_const_region singleton(const double& x) const;

	/*
	* region $(a_1, b_1]$ is considered "less than" $(a_2, b_2]$ if $a_1 < a_2$.
	*/
	bool operator<(const real_const_region& x) const;

	/*
	* region $(a_1, b_1]$ is considered "equal to" $(a_2, b_2]$ if $a_1 = a_2$
	* and $b_1 = b_2$.
	*/
	bool operator==(const real_const_region& x) const;

	/*
	* Set this region to be equal to `x`.
	*/
	const real_const_region& operator=(const real_const_region& x);

protected:
	double _a;
	double _b;
	dfdb _w;
	univariate_helper _helper;
	double _log_w_max;
	double _log_w_min;
	double _log_prob;
	optimizer _maxopt;
	optimizer _minopt;
	vws::midpoint _mid;
};

inline void real_const_region::set_w(const dfdb& w)
{
	_w = w;
}

inline void real_const_region::set_helper(const univariate_helper& helper)
{
	_helper = helper;
}

inline void real_const_region::set_maxopt(const optimizer& maxopt)
{
	_maxopt = maxopt;
}

inline void real_const_region::set_minopt(const optimizer& minopt)
{
	_minopt = minopt;
}

inline void real_const_region::set_mid(const vws::midpoint& mid)
{
	_mid = mid;
}

inline real_const_region::real_const_region(double a,
	const dfdb& w, const univariate_helper& helper,
	const optimizer& maxopt, const optimizer& minopt, const vws::midpoint& mid)
: _a(a), _b(a), _w(w), _helper(helper), _log_w_max(NAN), _log_w_min(NAN),
  _log_prob(NAN), _maxopt(maxopt), _minopt(minopt), _mid(mid)
{
	init();
}

inline real_const_region::real_const_region(double a, double b,
	const dfdb& w, const univariate_helper& helper,
	const optimizer& maxopt, const optimizer& minopt, const vws::midpoint& mid)
: _a(a), _b(b), _w(w), _helper(helper), _log_w_max(NAN), _log_w_min(NAN),
  _log_prob(NAN), _maxopt(maxopt), _minopt(minopt), _mid(mid)
{
	init();
}

inline void real_const_region::init()
{
	if (_a > _b) {
		// Invalid interval
		Rcpp::stop("a > b");
	} else if (_a < _b) {
		// Nontrivial interval (a,b].
		// Compute $P(a < X <= b)$ for $X \sim g$ on the log scale.
		_log_w_max = _maxopt(_w, _a, _b, true);
		_log_w_min = _minopt(_w, _a, _b, true);
		_log_prob = log_sub2_exp(_helper.p(_b, true, true), _helper.p(_a, true, true));
	} else {
		// Singleton interval (a,a].
		_log_w_max = _w(_a, true);
		_log_w_min = _w(_a, true);
		_log_prob = R_NegInf;
	}

	if ( (_log_w_max > 0) && std::isinf(_log_w_max) ) {
		Rcpp::stop("%s. %s",
			"Infinite maximum value found in optimize",
			"Cannot be used with real_const_region");
	}
}

inline double real_const_region::d_base(const double& x, bool log) const
{
	return _helper.d(x, log);
}

inline std::vector<double> real_const_region::r(unsigned int n) const
{
	// Generate a draw from $g_j$; i.e., the density $g$ truncated to this
	// region. Compute q((pb - pa) * u + pa) on the log scale.
	const Rcpp::NumericVector& u = Rcpp::runif(n);
	double log_pa = _helper.p(_a, true, true);
	const Rcpp::NumericVector& log_p = log_add2_exp(_log_prob + log(u), Rcpp::rep(log_pa, n));

	std::vector<double> out;
	for (unsigned int i = 0; i < n; i++) {
		out.push_back(_helper.q(log_p(i), true, true));
	}

	return out;
}

inline bool real_const_region::s(const double& x) const
{
	return _a < x && x <= _b;
}

inline double real_const_region::w(const double& x, bool log) const
{
	return _w(x, log);
}

inline double real_const_region::w_major(const double& x, bool log) const
{
	double out = _log_w_max;
	return log ? out : exp(out);
}

inline double real_const_region::w_minor(const double& x, bool log) const
{
	double out = _log_w_min;
	return log ? out : exp(out);
}

inline real_const_region real_const_region::merge(const real_const_region& x) const
{
	if (!is_mergeable(x)) {
		Rcpp::stop("Cannot merge these regions");
	}

	if (_a < x._a) {
		return real_const_region(_a, x._b, _w, _helper, _maxopt, _minopt, _mid);
	} else {
		return real_const_region(x._a, _b, _w, _helper, _maxopt, _minopt, _mid);
	}
}

inline double real_const_region::midpoint() const
{
	return _mid(_a, _b);
}

inline std::pair<real_const_region,real_const_region>
real_const_region::bifurcate() const
{
	return bifurcate(midpoint());
}

inline std::pair<real_const_region,real_const_region>
real_const_region::bifurcate(const double& x) const
{
	real_const_region r1(_a, x, _w, _helper, _maxopt, _minopt, _mid);
	real_const_region r2(x, _b, _w, _helper, _maxopt, _minopt, _mid);
	return std::make_pair(r1, r2);
}

inline real_const_region real_const_region::singleton(const double& x) const
{
	return real_const_region(x, _w, _helper, _maxopt, _minopt, _mid);
}

inline bool real_const_region::is_bifurcatable() const
{
	return true;
}

inline bool real_const_region::is_mergeable(const real_const_region& x) const
{
	return (this->_a == x._b || this->_b == x._a);
}


inline double real_const_region::xi_upper(bool log) const
{
	double out = _log_w_max + _log_prob;
	return log ? out : exp(out);
}

inline double real_const_region::xi_lower(bool log) const
{
	double out = _log_w_min + _log_prob;
	return log ? out : exp(out);
}

inline std::string real_const_region::description() const
{
	char buf[32];
	snprintf(buf, sizeof(buf), "(%g, %g]", _a, _b);
	return buf;
}

inline bool real_const_region::operator<(const real_const_region& x) const
{
	return _b <= x._a;
}

inline bool real_const_region::operator==(const real_const_region& x) const
{
	return _a == x._a && _b == x._b;
}

inline const real_const_region& real_const_region::operator=(const real_const_region& x)
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
	_mid = x._mid;
	return *this;
}

}

#endif

