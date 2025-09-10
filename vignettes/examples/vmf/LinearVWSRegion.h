#ifndef LINEAR_VWS_REGION_H
#define LINEAR_VWS_REGION_H

#include <Rcpp.h>
#include "vws.h"
#include "texp.h"
#include "target.h"
#include "fntl.h"

class LinearVWSRegion : public vws::Region<double>
{
private:
	double _a;
	double _b;
	double _kappa;
	double _d;
	double _beta0_min;
	double _beta1_min;
	double _beta0_max;
	double _beta1_max;

public:
	LinearVWSRegion(double a, double b, double kappa, double d);
	LinearVWSRegion(double a, double kappa, double d);
	void init();

	double midpoint() const;
	double d_base(const double& x, bool log = false) const;
	double w(const double& x, bool log = true) const;

	std::vector<double> r(unsigned int n) const;
	double d(const double& x, bool log = false) const;

	double get_lower() const { return _a; }
	double get_upper() const { return _b; }
	double get_d() const { return _d; }
	double get_kappa() const { return _kappa; }
	double get_beta0_min() const { return _beta0_min; }
	double get_beta1_min() const { return _beta1_min; }
	double get_beta0_max() const { return _beta0_max; }
	double get_beta1_max() const { return _beta1_max; }

	bool s(const double& x) const {
		return _a < x && x <= _b;
	}

	bool is_bifurcatable() const {
		double x = midpoint();
		return x > _a && x < _b;
	}

	double w_major(const double& x, bool log = true) const;
	double xi_upper(bool log = true) const;
	double xi_lower(bool log = true) const;

	std::pair<LinearVWSRegion,LinearVWSRegion> bifurcate() const {
		double x = midpoint();
		return bifurcate(x);
	}

	std::pair<LinearVWSRegion,LinearVWSRegion> bifurcate(const double& x) const;

	std::string description() const {
		char buf[32];
		sprintf(buf, "(%g, %g]", _a, _b);
		return buf;
	}

	LinearVWSRegion singleton(const double& x) const {
		return LinearVWSRegion(x, _kappa, _d);
	}

	void print() const {
		printf("Linear VWS Region (%g, %g]\n", _a, _b);
	}

	bool operator<(const LinearVWSRegion& x) const {
		return _a < x._a;
	}

	bool operator==(const LinearVWSRegion& x) const {
		return _a == x._a && _b == x._b;
	}

	const LinearVWSRegion& operator=(const LinearVWSRegion& x);
};

inline double LinearVWSRegion::w(const double& x, bool log) const
{
	double out = 0.5 * (_d - 3) * log1p(-std::pow(x, 2)) + std::log(-1 < x && x < 1);
	return log ? out : exp(out);
}

inline double LinearVWSRegion::d_base(const double& x, bool log) const
{
	return d_texp(x, _kappa, -1, 1, log);
}

inline double LinearVWSRegion::d(const double& x, bool log) const
{
	double kappa_j = _kappa + _beta1_max;
	return d_texp(x, kappa_j, _a, _b, log);
}

inline std::vector<double> LinearVWSRegion::r(unsigned int n) const
{
	double kappa_j = _kappa + _beta1_max;
	const auto& out = r_texp(n, kappa_j, _a, _b);
	return Rcpp::as<std::vector<double>>(out);
}

inline double LinearVWSRegion::w_major(const double& x, bool log) const
{
	double out = s(x) ? _beta0_max + _beta1_max * x : R_NegInf;
	return log ? out : exp(out);
}

inline double LinearVWSRegion::midpoint() const
{
	double out;

	if (std::isinf(_a) && std::isinf(_b) && _a < 0 && _b > 0) {
		// In this case, we have an interval (-Inf, Inf). Make a split at zero.
		out = 0;
	} else if (std::isinf(_a) && _a < 0) {
		// Left endpoint is -Inf. Split based on right endpoint.
		out = _b - std::fabs(_b) - 1;
	} else if (std::isinf(_b) && _b > 0) {
		// Right endpoint is Inf. Split based on left endpoint.
		out = _a + std::fabs(_a) + 1;
	} else {
		out = (_a + _b) / 2;
	}

	return out;
}

inline LinearVWSRegion::LinearVWSRegion(double a, double b, double kappa, double d)
: _a(a), _b(b), _kappa(kappa), _d(d), _beta0_min(), _beta1_min(),
  _beta0_max(), _beta1_max()
{
	init();
}

inline LinearVWSRegion::LinearVWSRegion(double a, double kappa, double d)
: _a(a), _b(a), _kappa(kappa), _d(d), _beta0_min(), _beta1_min(),
  _beta0_max(), _beta1_max()
{
	// Special handling for singleton sets.
	_beta0_min = a;
	_beta1_min = 0;
	_beta0_max = a;
	_beta1_max = 0;
}

inline void LinearVWSRegion::init()
{
	if (_a >= _b) {
		Rcpp::stop("a >= b: %g >= %g", _a, _b);
	}

	// First derivative of log-weight function
	const std::function<double(double)>& d_log_w =
	[&](double x) -> double {
		return -(_d - 3) * x / (1 - std::pow(x, 2));
	};

	const std::function<double(double)>& obj =
	[&](double x) -> double {
		double dx = d_log_w(x);
		return w(x, true) - x * dx + mgf_texp(dx, _kappa, _a, _b, true);
	};

	// Determine if log(w(x)) is convex, concave, or a constant. Majorization
	// and minorization choices will depend on this.
	if (_d > 3) {
		// log(w(x)) is concave. Any tangent line will lie above the curve.

		// For the majorizer, solve a minimization problem
		fntl::optimize_args args;
		args.fnscale = 1;
		const auto& opt_out = fntl::optimize_brent(obj, _a, _b, args);
		double c_star = opt_out.par;
		_beta0_max = w(c_star, true) - c_star*d_log_w(c_star);
		_beta1_max = d_log_w(c_star);

		// For the minorizer
       _beta1_min = (w(_b, true) - w(_a, true)) / (_b - _a);
       _beta0_min = w(_a, true) - _a*_beta1_min;
	} else if (_d < 3) {
		// log(w(x)) is convex. Make a line that passes through a and b.

		// For the majorizer
       _beta1_max = (w(_b, true) - w(_a, true)) / (_b - _a);
       _beta0_max = w(_a, true) - _a*_beta1_max;

		// For the minorizer, solve a maximization problem
		fntl::optimize_args args;
		args.fnscale = -1;
		const auto& opt_out = fntl::optimize_brent(obj, _a, _b, args);
		double c_star = opt_out.par;
		_beta0_min = w(c_star, true) - c_star*d_log_w(c_star);
		_beta1_min = d_log_w(c_star);
	} else {
		// log(w(x)) is constant with value zero
		_beta0_max = 0;
		_beta1_max = 0;
		_beta0_min = 0;
		_beta1_min = 0;
	}
}

inline double LinearVWSRegion::xi_upper(bool log) const
{
	double kappa_j = _kappa + _beta1_max;
	double lnc0 = n_texp(_kappa, -1, 1, true);
	double lncj = n_texp(kappa_j, _a, _b, true);
	double out = _beta0_max + lncj - lnc0;
	return log ? out : exp(out);
}

inline double LinearVWSRegion::xi_lower(bool log) const
{
	// Use the trivial minorizer here: integrate the original target density.
	double lnc0 = n_texp(_kappa, -1, 1, true);
	double lnc = ::n_target(_kappa, _d, true);
	double lp = ::integrate_target(_a, _b, _kappa, _d, true);
	double out = lp + lnc - lnc0;
	double log_xi_upper = xi_upper(true);
	if (log_xi_upper < out) {
		// This condition can happen numerically. If it occurs, just take lower
		// to be equal to upper.
		Rprintf("LinearVWSRegion: log_xi_lower (%g) <- log_xi_upper (%g)\n",
			out, log_xi_upper);
		out = log_xi_upper;
	}

	return log ? out : exp(out);
}

inline std::pair<LinearVWSRegion,LinearVWSRegion> LinearVWSRegion::bifurcate(const double& x) const
{
	LinearVWSRegion r1(_a, x, _kappa, _d);
	LinearVWSRegion r2(x, _b, _kappa, _d);
	return std::make_pair(r1, r2);
}

inline const LinearVWSRegion& LinearVWSRegion::operator=(const LinearVWSRegion& x)
{
	_a = x._a;
	_b = x._b;
	_kappa = x._kappa;
	_d = x._d;
	_beta0_min = x._beta0_min;
	_beta1_min = x._beta1_min;
	_beta0_max = x._beta0_max;
	_beta1_max = x._beta1_max;
	return *this;
}

#endif

